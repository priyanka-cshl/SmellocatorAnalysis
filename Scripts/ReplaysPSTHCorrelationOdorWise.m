MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat';
MySession = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/O3/O3_20211005_r0_processed.mat';

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

for i = 1:numel(MyUnits) % every unit
    subplot(6,1,i);  hold on
    MyCorrs = squeeze(MedianCorrs(i,:,:))';
    MySTDs = squeeze(STDCorrs(i,:,:))';
    bw1 = 0.2; bw2 = 0.22;
    
    for odor = 1:4
        whichones = find(OdorSequence==odor);
        xpts = (whichones*4)' - 3;
        if odor == 4
           bw1 = 0.9; bw2 = 1;
        end
            
        
        % OL-OL
        bar(xpts, MyCorrs(whichones,3),'Facecolor',Plot_Colors(['Odor',num2str(odor)]),...
            'BarWidth',bw1,'LineStyle','none');
        line(repmat(xpts,2,1), ...
            [(MyCorrs(whichones,3)' + MySTDs(whichones,3)'); (MyCorrs(whichones,3)' - MySTDs(whichones,3)')], ...
            'color', Plot_Colors(['Odor',num2str(odor)]), 'Linewidth', 2);
        
        %CL-OL
        bar(xpts, MyCorrs(whichones,1),'Edgecolor','k','Facecolor','none','BarWidth',bw2,'Linewidth', 2);
        line(repmat(xpts,2,1), ...
            [(MyCorrs(whichones,1) + MySTDs(whichones,1))'; (MyCorrs(whichones,1) - MySTDs(whichones,1))'], ...
            'color', [0 0 0],'Linewidth', 2);
        
        %PR-OL
        xpts = xpts + 1.5;
        plot(xpts, MyCorrs(whichones,4), 'o', ...
            'MarkerSize',6,'MarkerEdgecolor',[0.4 0.6 0.6],'MarkerFacecolor',[0.4 0.6 0.6]);
        line(repmat(xpts,2,1), ...
            [(MyCorrs(whichones,4) + MySTDs(whichones,4))'; (MyCorrs(whichones,4) - MySTDs(whichones,4))'], ...
            'color', [0.1 0.4 0.2], 'Linewidth', 2);

    end
end
            set(gca, 'XTick', xpts);
        xticklabels(num2str(1:21));
        xtickangle(gca,90);
%% Population PSTH
for i = 1:size(OpenLoopPSTH,1) % every replay
    popPSTH(:,i) = mean(squeeze(OpenLoopPSTH(i,:,:)),2);
    errPSTH(:,i) = std(squeeze(OpenLoopPSTH(i,:,:))');
end