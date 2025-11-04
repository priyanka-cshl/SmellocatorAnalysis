%% script to predict the sniff aligned spiking responses of given neuron 
%  post fitting with the poisson GLM

% sessionPath = 'T3_20250516_r0_processed.mat'; % path to the preprocessed session from PG
% GLMPaths = WhichComputerGLM;
% fullsessionPath = fullfile(GLMPaths.processedfile,sessionPath);
% load(sessionPath);
% 
% % sniff traces
% SniffTrace(:,1) = TracesOut.Timestamps{1};
% SniffTrace(:,2) = TracesOut.SniffsDigitized{1};

%% Step 1: Data Paths
myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

%% Get the sniff timestamps
%  Load Sniffs from the T3 session and make a vector with binarized sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
SniffTrace(:,1) = RespirationData(:,1); % Timestamps in 500hz base
SniffTrace(:,2) = 0;
% make a digitized trace
for i = 1:size(AllSniffs,1) % every sniff timestamp
    idx = AllSniffs(i,11:12);
    SniffTrace(idx(1):idx(2),2) = 1;
end

% downsample to 10ms binsize
PSTHbinsize = 10;
downsample  = PSTHbinsize/(1000/500);
TS_temp     = SniffTrace(:,1)';
TS_down     = TS_temp(1):PSTHbinsize/1000:TS_temp(end);
SniffTraceDS(:,1) = TS_down;
SniffTraceDS(:,2) = interp1(TS_temp,SniffTrace(:,2)',TS_down,'nearest');
SniffTraceDS(2:end,3) = diff(SniffTraceDS(:,2));
SniffTraceDS(SniffTraceDS(:,2)~=1,3) = 0;

%% Load the kernels
load('/mnt/data/sniffGLM/T3_20250516_r0_shared/SixteenOdorsFitOutput.mat');
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

%% for a given unit
% selectedUnits = [30 56];
whichUnit = 56;

myKernel    = kernels.pGLM(:,1,whichUnit);
myBaseline  = squeeze(baseline(1,1,whichUnit));
predictedFR = exp(conv(SniffTraceDS(:,3),myKernel,'full') + myBaseline);
predictedSpikes = PechePourPoisson(100*predictedFR,0.01);

%% Sniff timestamps
AllSniffs = SniffTraceDS(find(diff(SniffTraceDS(:,3))),1);
AllSniffs(1:end-1,2) = AllSniffs(2:end,1) - AllSniffs(1:end-1,1);

load(fullfile(myKsDir,'SessionDetails.mat'));
startTS = 0;
for i = 1:numel(Files.Samples)
    endTS = Files.Samples(i)/30000 + startTS; % in seconds
    whichsniffs = find( (AllSniffs(:,1)>=startTS) & (AllSniffs(:,1)<endTS) );
    AllSniffs(whichsniffs,3) = i; % session phase
    startTS = endTS;
end

AllSniffs = sortrows(AllSniffs,[3 2],'ascend');
thisUnitSpikes = SingleUnits(whichUnit).spikes;
[SpikeRaster, maxsniffs] = GetSniffLockedSpikes({AllSniffs}, thisUnitSpikes);
subplot(1,2,1)
plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);

[SpikeRaster, maxsniffs] = GetSniffLockedSpikes({AllSniffs}, predictedSpikes);
subplot(1,2,2)
plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);




