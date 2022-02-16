if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,'O3','O3_20211005_r0_processed.mat');

%% get the processed data loaded and assemble open loop traces
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

N = size(SingleUnits,2); % #units
MyUnits = (1:N);
%MyUnits = [58 35 34 55 21 8];

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

OdorSequence = OpenLoop.TTLs.OdorValve{1}(:,4);
if numel(OdorSequence)>size(TrialTS,1)
    OdorSequence(1,:) = [];
end

%% Split the long trace into odor-specific stretches 
% - only keep points from this trial's odorstart to next trial's odor start
for whichodor = 1:3
    whichones = find(OdorSequence==whichodor);
    MyIdx = [];
    for i = 1:numel(whichones)
        idx(1) = TrialTS(whichones(i),1); % odor start
        if whichones(i)<size(TrialTS,1)
            idx(2) = TrialTS(whichones(i)+1,1) - 1; % next trial odor start
        else
            idx(2) = TrialTS(whichones(i),3) + SampleRate; % 1 sec post trial off
        end
        MyIdx = horzcat(MyIdx,idx(1):idx(2));
    end
    [C{whichodor}] = ReplayCrossCorr2(OpenLoopPSTH(:,MyIdx,:),[1 10 5]);
    [Residuals{whichodor}] = ReplayResiduals(OpenLoopPSTH(:,MyIdx,:),[1 10 5]);
end
% also get correlation for the full trace
idx(1) = TrialTS(1,1);
idx(2) = TrialTS(end,3) + SampleRate;
[C{end+1}] = ReplayCrossCorr2(OpenLoopPSTH(:,idx(1):idx(2),:),[1 10 5]);
[Residuals{end+1}] = ReplayResiduals(OpenLoopPSTH(:,idx(1):idx(2),:),[1 10 5]);

%% Plotting stuff
for i = 1:N
    for j = 1:4
        withinOL(i,:,j) = C{j}(i).withinOL;
        withinPR(i,:,j) = C{j}(i).withinPR;
        CLvsOL(i,:,j) =  C{j}(i).CL2OL;
    end
end

%%
figure;
for j = 1:4
    subplot(3,4,j);
    histogram(withinOL(:,1,j),20,'Normalization','probability')
    subplot(3,4,j+4);
    histogram(withinPR(:,1,j),20,'Normalization','probability')
    subplot(3,4,j+8);
    histogram(CLvsOL(:,1,j),20,'Normalization','probability')
end