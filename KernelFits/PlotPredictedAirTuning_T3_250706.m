%% script to extract the sniff aligned spiking responses
%  of a given neuron 

%% Step 1: Data Paths
myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

%% Get the sniff timestamps
%  Load Sniffs from the T3 session and make a vector with binarized sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
load(fullfile(myKsDir,'SessionDetails.mat'));
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
Sniffs2Use{1}(11000:end,:) = [];
Sniffs2Use{1}(5000:8500,:) = [];

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
nRows = 6*2; nCols = 10; 

for thisunit = 1:numel(selectedUnits)
    
    if ~rem(thisunit-1,nCols);
        figure;
    end
    
    ts = (1:41)*0.05; % 20hz resolution
    ts_new = 0.02:0.02:ts(end);
    k = selectedUnits(thisunit)';
    KernelDS = interp1(ts,KernelsSmooth(k,:),ts_new);
    KernelDS(:,find(isnan(KernelDS))) = [];
    
    % convolve with the sniff trace
    PSTHout = conv(SniffTraceDS(:,2)',KernelDS,'full');
    PSTHout(PSTHout<0) = 0;
    %PSTHout(:,1:length(KernelDS)) = [];
    thisUnitSpikes = PechePourPoisson(1*PSTHout,PSTHbinsize/1000);

    % sort by duration only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},3,'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    whichplots = (1:10:(nRows*nCols*0.5)) + rem(thisunit-1,nCols);
    subplot(nRows+1,nCols,whichplots);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    
    % sort by session phase and duration only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},[8 3],'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    whichplots = whichplots + (nRows*nCols*0.5);
    subplot(nRows+1,nCols,whichplots);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);

    subplot(nRows+1,nCols,whichplots(end)+10);
    plot(KernelDS);

    if ~rem(thisunit,nCols) || thisunit == numel(selectedUnits)
        set(gcf,'Position',[2015 535 1748 323]);
        if savefigs
            set(gcf,'Color','w');
            set(gcf,'renderer','Painters');
            print(['/mnt/data/simulation/SelectedPredictedAirResponse',num2str(thisunit),'.eps'],'-depsc','-tiff','-r300','-painters');
            close(gcf);
        end
    end
    
    
end
