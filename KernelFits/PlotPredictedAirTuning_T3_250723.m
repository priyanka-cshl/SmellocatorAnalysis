%% script to extract the sniff aligned spiking responses
%  of a given neuron 

%% Step 1: Data Paths
myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

%% Get the sniff timestamps
%  Load Sniffs from the T3 session and make a vector with binarized sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
load(fullfile(myKsDir,'SessionDetails.mat'));

if exist('CuratedSniffTimestamps','var')
    AllSniffs = CuratedSniffTimestamps(find(CuratedSniffTimestamps(:,8)>0),:);
    AllSniffs(:,11:12) = AllSniffs(:,8:9);
    AllSniffs(:,4:7) = 0;
    AllSniffs(:,3) = AllSniffs(:,3) - AllSniffs(:,1);
end

% add odor TTL info if available
if exist(fullfile(myKsDir,'quickprocessOdorTTLs.mat'))
    load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings');
    
    % add the column infos to the sniffs
    for i = 1:size(TTLs.Trial,1) % everty trial
        t = TTLs.Trial(i,7:10); % odor1 start, stop, purge start, stop
        t(end+1) = t(end) + TTLs.Trial(i,9)-TTLs.Trial(i,8); % add a post-purge period the same as the first pulse
        TS = [t(1:end-1)' t(2:end)' [1 3 2 4]'];
        for j = 1:size(TS,1)
            whichsniffs = find( (AllSniffs(:,1)>=TS(j,1)) & (AllSniffs(:,1)<TS(j,2)) );
            AllSniffs(whichsniffs,5) = TTLs.Trial(i,4); % odor identity
            AllSniffs(whichsniffs,6) = TS(j,3);
        end
    end
    
end

startTS = 0;
for i = 1:numel(Files.Samples)
    endTS = Files.Samples(i)/30000 + startTS; % in seconds
    whichsniffs = find( (AllSniffs(:,1)>=startTS) & (AllSniffs(:,1)<endTS) );
    AllSniffs(whichsniffs,8) = i; % session phase
    startTS = endTS;
end

% group sniffs by odors
% there are no air OFF sniffs here
AllSniffs(:,4) = 1; % air is always on
[ParsedSniffs, StimulusList] = ParseSniffsByStimuli(AllSniffs, 'SortBy', 3);
Sniffs2Use{1} = ParsedSniffs{2}; % Air sniffs
% Sniffs2Use{1}(:,8) = 1;
% Sniffs2Use{1}(11000:end,:) = [];
% Sniffs2Use{1}(5000:8500,:) = [];

%% SniffTrace for convolving
SniffTrace(:,1) = RespirationData(:,1); % Timestamps in 500hz base
SniffTrace(:,2) = 0;

for i = 1:size(AllSniffs,1) % every sniff timestamp
    idx = AllSniffs(i,11:12);
    SniffTrace(idx(1):idx(2),2) = 1;
end

PSTHbinsize = 20; 
% downsample to 20ms binsize
downsample  = PSTHbinsize/(1000/500);
TS_temp     = SniffTrace(:,1)';
TS_down     = TS_temp(1):PSTHbinsize/1000:TS_temp(end);

SniffTraceDS(:,1) = TS_down;
SniffTraceDS(:,2) = interp1(TS_temp,SniffTrace(:,2)',TS_down,'nearest');

%% load MT kernels
load('/mnt/data/simulation/MTKernels.mat');
selectedUnits = [5 6 14 15 33 160 194 195 196 201]; % 134 167

%%
savefigs = 1;
ntypes = 3;
nRows = 6*ntypes; nCols = 10; 

for thisunit = 1:numel(selectedUnits) %size(KernelsSmooth,1) %
    
    if ~rem(thisunit-1,nCols);
        figure;
    end
    baseline = 5;

    ts = (1:41)*0.025; % 20hz resolution
    ts_new = 0.02:0.02:ts(end);
    k = selectedUnits(thisunit)';
%     k = thisunit;
    KernelDS = interp1(ts,KernelsSmooth(k,:),ts_new);
    KernelDS(:,find(isnan(KernelDS))) = [];
    KernelDS = KernelDS/3;

    % convolve with the sniff trace
    PSTHMain = conv(SniffTraceDS(:,2)',KernelDS,'full');

    PSTHout = PSTHMain + baseline ;
    PSTHout(PSTHout<0) = 0;
    %PSTHout(:,1:length(KernelDS)) = [];
    thisUnitSpikes = PechePourPoisson(1*PSTHout,PSTHbinsize/1000);

    % sort by duration only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},3,'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    whichplots = (1:10:(nRows*nCols/ntypes)) + rem(thisunit-1,nCols);
    subplot(nRows+1,nCols,whichplots);
    %plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    h1 = plot(nan,nan,'.k','Markersize', 0.5);
    h1.XData = SpikeRaster{1}(:,1)'; 
    h1.YData = SpikeRaster{1}(:,2)';
    hold on
    plot(Sniffs2Use{1}(:,1)*0,1:size(Sniffs2Use{1},1),'r');
    plot(Sniffs2Use{1}(:,3),1:size(Sniffs2Use{1},1),'r');
    
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},[8 3],'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    whichplots = whichplots + (nRows*nCols/ntypes);
    subplot(nRows+1,nCols,whichplots);
    %plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    h2 = plot(nan,nan,'.k','Markersize', 0.5);
    h2.XData = SpikeRaster{1}(:,1)';
    h2.YData = SpikeRaster{1}(:,2)';
    hold on
    plot(Sniffs2Use{1}(:,1)*0,1:size(Sniffs2Use{1},1),'r');
    plot(Sniffs2Use{1}(:,3),1:size(Sniffs2Use{1},1),'r');

    Sniffs2Use{1} = sortrows(Sniffs2Use{1},[8 1],'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);

    whichplots = whichplots + (nRows*nCols/3);
    subplot(nRows+1,nCols,whichplots);
    %plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    h3 = plot(nan,nan,'.k','Markersize', 0.5);
    h3.XData = SpikeRaster{1}(:,1)';
    h3.YData = SpikeRaster{1}(:,2)';

    subplot(nRows+1,nCols,whichplots(end)+10);
    plot(KernelDS);

    if ~rem(thisunit,nCols) || thisunit == numel(selectedUnits)
        set(gcf,'Position',[2015 6 1748 989]);
        if savefigs
            set(gcf,'Color','w');
            set(gcf,'renderer','Painters');
            print(['/mnt/data/simulation/HalfFastSelectedPredictedAirResponse',num2str(thisunit),'.eps'],'-depsc','-tiff','-r300','-painters');
            close(gcf);
        end
    end
    
    
end

%%
% figure; 
% nToT = 110;
% for thisunit = 1:size(KernelsSmooth,1) %numel(selectedUnits)
%     
%     if ~rem(thisunit-1,nToT);
%         figure;
%     end
%     
%     ts = (1:41)*0.025; % 20hz resolution
%     ts_new = 0.02:0.02:ts(end);
%     %k = selectedUnits(thisunit)';
%     k = thisunit;
%     KernelDS = interp1(ts,KernelsSmooth(k,:),ts_new);
%     KernelDS(:,find(isnan(KernelDS))) = [];
% 
%     subplot(11,10,k);
%     plot(KernelDS);
% end

