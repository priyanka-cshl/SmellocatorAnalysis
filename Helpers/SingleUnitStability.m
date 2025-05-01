function [cluster] = SingleUnitStability(myKsDir, varargin)
% function to load spike waveforms for each unit an assess unit stability over time

%% extract inputs
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotwaveforms', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('savefigures', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotWaveforms = params.Results.plotwaveforms;
saveFigures = params.Results.savefigures;
plotDrift = 1;

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

if saveFigures
    if ~exist(fullfile(myKsDir,'ClusterMaps'),'dir')
        mkdir(fullfile(myKsDir,'ClusterMaps'));
    end
end

for mycluster = 1:length(sp.cids) % for each cluster

    disp(mycluster);

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
    [waveForms{1}, meanWFs{1}] = extractWaveforms(whichIndices, mmf); % numWF x nChannels x WinSize
    
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
        [waveForms{seg+1}, meanWFs{seg+1}] = extractWaveforms(whichIndices, mmf); % numWF x nChannels x WinSize
    end

    cluster(mycluster).WF = waveForms;
    cluster(mycluster).meanWF = meanWFs;

    if plotWaveforms || saveFigures

        % plot waveform map
        figure;
        nTetrodePairs = ceil(nChInMap/8);
        subplotRows = nTetrodePairs + 4;
        % determine YLim by taking the waveform span on the strongest channel
        myLims = [ min(min(meanWFs{1}(otherChannels+1,:,1))) ...
            max(max(meanWFs{1}(otherChannels+1,:,1))) ];
        myLims = 100*[floor(myLims(1)/100) ceil(myLims(2)/100)];

        % plotting
        for whichChannel = 1:nChInMap
            whichPlot = whichChannel;
            subplot(subplotRows,8,whichPlot);
            meanWF = squeeze(meanWFs{1}(whichChannel,:,1));
            stdWF = squeeze(meanWFs{1}(whichChannel,:,2));
            if ismember(whichChannel-1,otherChannels)
                MyShadedErrorBar([],meanWF,stdWF,'r',{},0.5);
            else
                MyShadedErrorBar([],meanWF,stdWF,'b',{},0.5);
            end

            if plotDrift
            % add the three recording segments
            offset = diff(myLims)/4;
            for seg = 1:3
                meanWF = squeeze(meanWFs{seg+1}(whichChannel,:,1)) - (offset*seg);
                stdWF = squeeze(meanWFs{seg+1}(whichChannel,:,2));
                if ismember(whichChannel-1,otherChannels)
                    MyShadedErrorBar([],meanWF,stdWF,'r',{},0.5);
                else
                    MyShadedErrorBar([],meanWF,stdWF,'b',{},0.5);
                end
            end
            end

            set(gca,'Box','off','Color','none','XColor','none','YColor','none',...
                'YLim',[myLims(1)-offset*3 myLims(2)],'XLim',[0 diff(wfWin)],'XTick',[],'YTick',[]);

            set(gca,'Box','off','Color','none','XColor','none','YColor','none',...
                'YLim',myLims,'XLim',[0 diff(wfWin)],'XTick',[],'YTick',[]);
            
        end

        % add the ISI and corr plots
        % map out the plotting grid
        foo = reshape(1:subplotRows*8,8,subplotRows)';
        ISI_subplots  = foo((nTetrodePairs+1):(end-1),1:4);
        Corr_subplots = foo((nTetrodePairs+1):(end-1),5:8);
        Amp_plots = foo(end,1:8);

        % plot the ISI histogram
        ISIs = 1000*diff(allspikes); % in milliseconds
        subplot(subplotRows,8,ISI_subplots(:));
        H = histogram(ISIs,[0:1:50 Inf]);
        myHist = H.Values(1:end-1);
        maxY = 100*ceil(max(myHist)/100);
        cla
        rectangle('Position',[0 0 2 maxY],'EdgeColor', 'none', 'FaceColor', 0.7*[1 1 1]);
        hold on
        bar(myHist,1,'EdgeColor','none');
        set(gca,'Color','none','YLim',[0 maxY],'XTick',[],'YTick',[]);
        title(['Cluster# ',num2str(cluster(mycluster).id)])

        
        % plot the autocorrelation
        SpikeTrain = zeros(ceil(max(1000*allspikes)),1);
        SpikeTrain(round(1000*allspikes),1) = 1;
        [r,lags] = xcorr(SpikeTrain,50,'coeff');
        r(lags==0) = NaN;
        subplot(subplotRows,8,Corr_subplots(:));
        bar(r,1,'EdgeColor','none');
        set(gca,'Color','none','XTick',[],'YTick',[]);
        title(['Amplitudes ',mat2str(myLims)])

        % plot the spike amplitudes over time
        subplot(subplotRows,8,Amp_plots(:));
        plot(cluster(mycluster).spikes, cluster(mycluster).Amplitudes,'.','MarkerSize', 4);
        axLims = get(gca,'YLim')/10;
        axLims(1) = floor(axLims(1))*10;
        axLims(2) = ceil(axLims(2))*10;
        set(gca,'Color','none','YLim',axLims,'YTick',axLims,'XLim',[0 ceil(recDuration)],'TickDir','out');
        
        if saveFigures
            % print to pdf
            set(gcf,'Color','w');
            set(gcf,'renderer','Painters');
            if mycluster == 1
                if nTetrodePairs == 8
                    figPosition = [0.7, 0.2, 0.2, 0.9];
                else
                    figPosition = [0.7, 0.2, 0.2, 0.7];
                end            
                set(gcf, 'Units', 'Normalized', 'OuterPosition', figPosition);
                exportgraphics(gcf, ...
                    fullfile(myKsDir,'ClusterMaps','UnitSummary.pdf'),...
                    'ContentType','vector');
            else
                set(gcf, 'Units', 'Normalized', 'OuterPosition', figPosition);
                exportgraphics(gcf, ...
                    fullfile(myKsDir,'ClusterMaps','UnitSummary.pdf'),...
                    'ContentType','vector','Append',true);
            end
            close;
        end
    end
    
end

save(fullfile(myKsDir,'ClustersFull.mat'),"cluster",'-v7.3');


    function [wfOut, meanwfOut] = extractWaveforms(mySpikeIndices, mmf)
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
            meanwfOut(:,:,1) = squeeze(mean(wfOut,1,'omitnan')); % average across all spikes
            if size(wfOut,1) > 1
                meanwfOut(:,:,2) = squeeze(std(double(wfOut),1,'omitnan'));
            end
    end
end