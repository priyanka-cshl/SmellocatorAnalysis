function [cluster] = SingleUnitStability(myKsDir, varargin)
% function to load spike waveforms for each unit an assess unit stability over time

%% extract inputs
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotwaveforms', true, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotWaveforms = params.Results.plotwaveforms;

%% add repos if needed
Paths = WhichComputer();
addpath(genpath([Paths.Code,filesep,'open-ephys-analysis-tools']));
addpath(genpath([Paths.Code,filesep,'afterphy']));
addpath(genpath([Paths.Code,filesep,'spikes']));
addpath(genpath([Paths.Code,filesep,'npy-matlab']));
addpath(genpath([Paths.Code,filesep,'MatlabUtils']));

%% Load data from kilosort/phy
sp = loadKSdirPriyanka(myKsDir);
% sp.st are spike times in seconds (for all spikes)
% sp.clu are cluster identities (for all spikes)
% sp.cids is list of unqiue clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted

%% for waveform extraction
fileName = fullfile(myKsDir,sp.dat_path); % binary .dat file that contains the raw data
filenamestruct = dir(fileName);
nChannels = sp.n_channels_dat;
if exist(fullfile(myKsDir,"SessionDetails.mat"))==2
    temp = load(fullfile(myKsDir,"SessionDetails.mat"));
    nSamples = temp.Files.Samples;
else
    dataTypeNBytes = numel(typecast(cast(0, sp.dtype), 'uint8')); % determine number of bytes per sample
    nSamples = filenamestruct.bytes/(nChannels*dataTypeNBytes);  % Number of samples per channel
end
mmf = memmapfile(fileName, 'Format', {sp.dtype, [nChannels nSamples], 'x'});
chMap = readNPY(fullfile(myKsDir, 'channel_map.npy'))+1; % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

% initializations for reading waveforms
wfWin = [-40 41];          % Number of samples before and after spiketime to include in waveform
nWf = 2000;                % Maximum Number of waveforms per unit to pull out
wfNSamples = length(wfWin(1):wfWin(end));

recDuration = max(sp.st);
nSegments   = 3;
segDuration = recDuration/nSegments;
timeSegments = repmat([0 segDuration],nSegments,1) + ...
                segDuration*((1:nSegments)-1)';

%% get all clusters and reorder by primary channel
[~,clusterorder] = sort(sp.channels);

for mycluster = 1:length(sp.cids) % for each cluster

    whichcluster = clusterorder(mycluster);
    
    % get all spiketimes (in seconds)
    allspikes = sp.st(sp.clu==sp.cids(whichcluster));
    primeChannel = sp.channels(whichcluster);
    otherChannels = floor(primeChannel/8)*8 + (1:8) -1;

    
    % which tetrode
    tetrode = floor(primeChannel/4)+1 + rem(primeChannel,4)/10;
    
    % Outputs
    cluster(mycluster).id = sp.cids(whichcluster);
    cluster(mycluster).tetrode = tetrode;
    cluster(mycluster).spikes = allspikes;
    cluster(mycluster).spikecount = numel(allspikes);
    cluster(mycluster).quality = sp.cgs(whichcluster);
    [fpRate, numViolations] = ISIViolations(allspikes, 1/30000, 0.002);
    cluster(mycluster).ISIquality = [round(fpRate,2,'significant'), round(numViolations/(numel(allspikes)-1),2,'significant')];
    cluster(mycluster).Amplitudes = sp.tempScalingAmps(sp.clu==sp.cids(whichcluster));
    
    % prep for waveform extraction
    % exclude any spikes earlier than the window limit 
    allspikes((allspikes<=abs(wfWin(1)/sp.sample_rate)),:) = [];
    allspikes((allspikes>=(recDuration-wfWin(2)/sp.sample_rate)),:) = [];
    % adjust number of waveforms to pull out as needed

    % for the whole session
    numWF       = min(nWf, numel(allspikes));
    permSpikes  = randperm(numel(allspikes));
    whichSpikes = allspikes(permSpikes(1:numWF));
    % convert to samples
    whichIndices = whichSpikes*sp.sample_rate + wfWin;
    whichIndices = sortrows(whichIndices,1);
    waveForms{1} = extractWaveforms(whichIndices, mmf); % numWF x nChannels x WinSize
    
    % for each time segment 
    % get total spike count in each segment
    for seg = 1:size(timeSegments,1)
        thisSegmentSpikes = allspikes(find(allspikes>timeSegments(seg,1) & allspikes<=timeSegments(seg,2)));
        cluster(mycluster).spikecount(seg+1) = numel(thisSegmentSpikes);
        numWF       = min(nWf, numel(thisSegmentSpikes));
        permSpikes  = randperm(numel(thisSegmentSpikes));
        whichSpikes = thisSegmentSpikes(permSpikes(1:numWF));
        % convert to samples
        whichIndices = whichSpikes*sp.sample_rate + wfWin;
        whichIndices = sortrows(whichIndices,1);
        waveForms{seg+1} = extractWaveforms(whichIndices, mmf); % numWF x nChannels x WinSize
    end
    
    if savefigs
        
        if ~exist(fullfile(myKsDir,'ClusterMaps'),'dir')
            mkdir(fullfile(myKsDir,'ClusterMaps'));
        end
        
        % plot waveform map
        figure;
        nTetrodes = 8 + 8*(sp.n_channels_dat > 40);
        % determine YLim by taking the waveform span on the strongest channel
        myLims = [min(min(wf.waveFormsMean(1,(tetrode-1)*4+[1 2 3 4],:))) ...
            max(max(wf.waveFormsMean(1,(tetrode-1)*4+[1 2 3 4],:))) ];
        myLims = 100*[floor(myLims(1)/100) ceil(myLims(2)/100)];
        for whichChannel = 1:sp.n_channels_dat
            whichPlot = 8*(sp.ycoords(whichChannel)/1000 - 1) + sp.xcoords(whichChannel);
            subplot(nTetrodes,8,whichPlot);
            %plot(squeeze(wf.waveFormsMean(1,whichChannel,:)));
            meanWF = squeeze(wf.waveFormsMean(1,whichChannel,:));
            stdWF = nanstd(squeeze(wf.waveForms(1,:,whichChannel,:)),1)';
            if ceil(whichChannel/4)==tetrode
                switch sp.cgs(whichcluster)
                    case 2 % good unit
                        MyShadedErrorBar([],meanWF,stdWF,'r');
                    otherwise
                        MyShadedErrorBar([],meanWF,stdWF,'b');
                end
                
            else
                MyShadedErrorBar([],meanWF,stdWF,'k');
            end
            set(gca,'Box','off','Color','none','XColor','none','YColor','none',...
                'YLim',myLims,'XLim',[0 diff(gwfparams.wfWin)],'XTick',[],'YTick',[]);
        end
        
        % plotting Grid
        foo = reshape(1:nTetrodes*8,8,nTetrodes)';
        if nTetrodes <= 8
            WF_subplots   = foo(1:2,5:8);
            ISI_subplots  = foo(3:5,5:8);
            Corr_subplots = foo(6:8,5:8);
        else
            WF_subplots   = foo(1:4,5:8);
            ISI_subplots  = foo(5:10,5:8);
            Corr_subplots = foo(11:16,5:8);
        end
        
        % plot individual waveforms
        myLims = [];
        for n = 1:4
            whichWire = 4*(tetrode-1) + n;
            subplot(nTetrodes,8,WF_subplots(:,n));
            plot(squeeze(wf.waveForms(1,:,whichWire,:))','r');
            myLims = vertcat(myLims,get(gca,'YLim'));
        end
        for n = 1:4
            subplot(nTetrodes,8,WF_subplots(:,n));
            set(gca,'YLim',[2*min(myLims(:,1)) 2*max(myLims(:,2))],...
                'Box','off','Color','none','XColor','none','YColor','none',...
                'XLim',[0 diff(gwfparams.wfWin)],'XTick',[],'YTick',[]);
        end
        
        % plot the ISI histogram
        ISIs = 1000*diff(allspikes); % in milliseconds
        subplot(nTetrodes,8,ISI_subplots(:));
        H = histogram(ISIs,[0:1:50 Inf]);
        myHist = H.Values(1:end-1);
        bar(myHist,1,'EdgeColor','none');
        set(gca,'Color','none','YLim',[0 100*ceil(max(myHist)/100)],'XTick',[],'YTick',[]);
        title(['Cluster# ',num2str(cluster(mycluster).id)])
        
        % plot the autocorrelation
        SpikeTrain = zeros(ceil(max(1000*allspikes)),1);
        SpikeTrain(round(1000*allspikes),1) = 1;
        [r,lags] = xcorr(SpikeTrain,50,'coeff');
        r(lags==0) = NaN;
        subplot(nTetrodes,8,Corr_subplots(:));
        bar(r,1,'EdgeColor','none');
        set(gca,'Color','none','XTick',[],'YTick',[]);
        
        % print to pdf
        set(gcf,'renderer','Painters');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.9]);
        print(fullfile(myKsDir,'ClusterMaps',['Cluster',num2str(cluster(mycluster).id),'.eps']),...
            '-depsc','-tiff','-r300','-painters');
        close;
    end
    
end

% sort by tetrodes
[~,sortorder] = sort(arrayfun(@(x) x.tetrode, cluster));
cluster = cluster(sortorder);

if ~nofuss
    cluster(1).spikescaling = sp.tempScalingAmps;
    cluster(1).clusterscalingorder = sp.clu;

    disp(['found ',num2str(mycluster),' units']);
end

    function [wfOut] = extractWaveforms(mySpikeIndices, mmf)
        % convert to samples
        for thisSpike = 1:size(mySpikeIndices,1)
            % every spike
            curSpikeIndex = int32(mySpikeIndices(thisSpike,:));
%             disp(curSpikeIndex);
%             if ~isinteger(curSpikeIndex(2))
%                 keyboard;
%             end
            wfOut(thisSpike,:,:) = mmf.Data.x(:,curSpikeIndex(1):curSpikeIndex(2));
        end
    end
end