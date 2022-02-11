MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat';

%% get the processed data loaded and assemble open loop traces
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

N = size(SingleUnits,2); % #units
MyUnits = (1:N);
MyUnits = [58 35 34 55 21 8];

[OpenLoopTraces,OpenLoopTimestamps,OpenLoopPSTH,~] = ...
    ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, 'whichunits', MyUnits);

% OpenLoopTraces [:,7,1+#activereplays]
% cols - lever, motor, sniffs, licks, rewards, trial+odor(-ve), TargetZone

%% Parse continuous traces into trials
Trial = OpenLoopTraces(:,6);
Trial(Trial<0) = 0;
Odor = abs(OpenLoopTraces(:,6));
TrialTS =  horzcat( find(diff(Odor)>0), find(diff(Trial)>0), find(diff(Trial)<0)); % Odor ON, Trial ON, Trial OFF
TrialTS(:,4) = Odor(TrialTS(:,2)); % which odor
TrialTS(:,5) = OpenLoopTraces(TrialTS(:,2),7); % which Target Zone

for t = 1:size(TrialTS,1) % every subtrial
    idx(1) = TrialTS(t,1); % odor start
    if t < size(TrialTS,1)
        idx(2) = TrialTS(t+1,1) - 1; % next trial odor start
    else
        idx(2) = TrialTS(t,3) + SampleRate; % 1 sec post trial off
    end
    [C{t}, ~] = ReplayCrossCorr(OpenLoopPSTH(:,idx(1):idx(2),:),[1 10 5]);
end
% also get correlation for the full trace
[C{t+1}, g] = ReplayCrossCorr(OpenLoopPSTH(:,TrialTS(1,1):idx(2),:),[1 10 5]);

%% Population PSTH
for i = 1:size(OpenLoopPSTH,1) % every replay
    popPSTH(:,i) = mean(squeeze(OpenLoopPSTH(i,:,:)),2);
    errPSTH(:,i) = std(squeeze(OpenLoopPSTH(i,:,:))');
end