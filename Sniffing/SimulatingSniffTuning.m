%% Settings
PSTHbinsize = 20; 

%% Load Sniffs from the T3 session and make a vector with binarized sniffs
load('/Users/Priyanka/Desktop/LABWORK_II/Data/SimulatingSniffTuning/quickprocesssniffs.mat');
SniffTrace(:,1) = RespirationData(:,1); % Timestamps in 500hz base
SniffTrace(:,2) = 0;

for i = 1:size(AllSniffs,1) % every sniff timestamp
    idx = AllSniffs(i,11:12);
    SniffTrace(idx(1):idx(2),2) = 1;
end

% figure; 
% plot(RespirationData(:,1),RespirationData(:,3));
% hold on
% plot(SniffTrace(:,1),SniffTrace(:,2));

% downsample to 20ms binsize
downsample  = PSTHbinsize/(1000/500);
TS_temp     = SniffTrace(:,1)';
TS_down     = TS_temp(1):PSTHbinsize/1000:TS_temp(end);

SniffTraceDS(:,1) = TS_down;
SniffTraceDS(:,2) = interp1(TS_temp,SniffTrace(:,2)',TS_down,'nearest');

%% Load Kernels
load('/Users/Priyanka/Desktop/LABWORK_II/Data/SimulatingSniffTuning/MTKernels.mat');

for k = 1:size(KernelsSmooth,1) % every kernel
    ts = (1:41)*0.05; % 20hz resolution
    ts_new = 0.02:0.02:ts(end);
    KernelDS = interp1(ts,KernelsSmooth(i,:),ts_new);
    KernelDS(:,find(isnan(KernelDS))) = [];

    % convolve with the sniff trace
    PSTHout = conv(SniffTraceDS(:,2)',KernelDS,'same');
    PSTHout(PSTHout<0) = 0;
    thisUnitSpikes = PechePourPoisson(PSTHout,0.020);
    SelectedSniffs{1} = AllSniffs;
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(SelectedSniffs, thisUnitSpikes);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);

end


