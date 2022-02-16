MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat';
%MySession = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/O3/O3_20211005_r0_processed.mat';

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

%% Plotting - distributions - looks ugly
for i = 1:numel(MyUnits) % every unit
    figure;
    for x = 1:21 % each trial separately and then whole session
        subplot(4,6,x);  hold on
        % get histograms
        h1 = histogram(C{x}(find(strcmp(g,'OL-OL')),i),[-1:0.1:1],...
            'Normalization','probability', ...
            'EdgeColor','none','FaceColor',Plot_Colors('Paletton',[2 2]));
        h2 = histogram(C{x}(find(strcmp(g,'PR-PR')),i),[-1:0.1:1],...
            'Normalization','probability', ...
            'EdgeColor','none','FaceColor',Plot_Colors('Paletton',[3 2]));
        h3 = histogram(C{x}(find(strcmp(g,'CL-OL')),i),[-1:0.1:1],...
            'Normalization','probability', ...
            'EdgeColor','k','FaceColor',Plot_Colors('Paletton',[2 1]));
        h4 = histogram(C{x}(find(strcmp(g,'CL-PR')),i),[-1:0.1:1],...
            'Normalization','probability', ...
            'EdgeColor','k','FaceColor',Plot_Colors('Paletton',[3 1])); 
    end
end

%% Plotting - just means and errors for each trial
U = unique(g);
for x = 1:size(TrialTS,1)+1 % each trial separately and then whole session
    for c = 1:length(U)
        MedianCorrs(:,c,x) = median(C{x}(find(strcmp(g,U{c})),:))';
        STDCorrs(:,c,x)    = std(C{x}(find(strcmp(g,U{c})),:))';
    end
end

OdorSequence = OpenLoop.TTLs.OdorValve{1}(:,4);
if numel(OdorSequence)>size(TrialTS,1)
    OdorSequence(1,:) = [];
end
OdorSequence(end+1) = 4;
%%

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