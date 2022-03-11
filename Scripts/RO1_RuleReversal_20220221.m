%% for PCX3
SessionPath = 'PCX3/PCX3_20210602_r0_processed.mat';

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'TrialInfo', 'SampleRate', 'startoffset', 'TargetZones');

% find trial indices where the block switches
BlockSwitches(:,2) = [find(diff(TrialInfo.TransferFunctionLeft)); TrialInfo.TrialID(end)];
BlockSwitches(:,1) = [1; BlockSwitches(1:end-1,2)+1];

nBlocks = size(BlockSwitches,1);
figure;
%% behavior plot
for i = 1:nBlocks-1
    subplot(nBlocks,1,i);
    whichTrials = [BlockSwitches(i,2) + [-9:1:0]'; BlockSwitches(i+1,1) + [0:1:9]'];
    for j = 1:numel(whichTrials)
        Traces.TargetZone(whichTrials(j)) = {TrialInfo.TargetZoneType(whichTrials(j)) + ...
                0*Traces.Trial{whichTrials(j)}};
        if TrialInfo.TransferFunctionLeft(whichTrials(j))
            Traces.TargetZone(whichTrials(j)) = {-Traces.TargetZone{whichTrials(j)}};
        end
    end
    
    [TracesOut] = ConcatenateTraces(Traces, whichTrials, SampleRate*startoffset);
    timestamps = (1:size(TracesOut.Lever{1},1))'/SampleRate;
    foo = TracesOut.TargetZone{1};
    TZ = [TargetZones(abs(foo),1) TargetZones(abs(foo),3)];
    TZ(foo<0,:) = -TZ(foo<0,:);
    Trial = TracesOut.Trial{1};
    
    PlotBehavior(timestamps,TracesOut.Lever{1},TracesOut.Sniffs{1},...
    TracesOut.Licks{1},TracesOut.Rewards{1},Trial,TZ, 5);
end


%% Success rate plot
binsize = 10; SuccessRate = [];
for i = binsize:binsize:TrialInfo.TrialID(end)
    x1 = i - binsize + 1;
    x2 = i;
    SuccessRate(i/binsize,1) = mean(TrialInfo.Success(x1:x2));
end
subplot(nBlocks,1,nBlocks); hold on
plot(SuccessRate,'o-');


SessionPath = 'PCX3/PCX3_20210610_r0_processed.mat';

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'TrialInfo', 'SampleRate', 'startoffset', 'TargetZones');
binsize = 10; SuccessRate = [];
for i = binsize:binsize:TrialInfo.TrialID(end)
    x1 = i - binsize + 1;
    x2 = i;
    SuccessRate(i/binsize,1) = mean(TrialInfo.Success(x1:x2));
end
plot(SuccessRate,'o-r');