%% paths
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,'O3','O3_20211005_r0_processed.mat');
%MySession = fullfile(datapath,'O8','O8_20220702_r0_processed.mat');
%MySession = fullfile(datapath,'Q4','Q4_20221112_r0_processed.mat');
%MySession = fullfile(datapath,'Q9','Q9_20221119_r0_processed.mat');

%% get the processed data loaded and assemble open loop traces
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

N = size(SingleUnits,2);
MyUnits = (1:N);

%% sort units by tetrodes and get open loop rasters and PSTHs
foo = cell2mat(arrayfun(@(x) [x.tetrode; x.id], SingleUnits, 'UniformOutput', false))';
[~, MyUnits] = sortrows(foo,[1 2]);

[OpenLoopTraces,OpenLoopTimestamps,OpenLoopPSTH,OpenLoopRaster] = ...
        ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'whichunits', MyUnits, 'PSTHsmooth', 100, ...
        'plotfigures', 0);

num_OL = numel(find(ReplayTTLs.TrialID<=TrialInfo.TrialID(end))); % active replays
num_PR = numel(ReplayTTLs.TrialID) - num_OL; % passive replays
reps_per_condition = [1 num_OL num_PR]; % CL, OL, PR

%% using the whole PSTH or just select timepoints eg. specific odor for comparison across conditions

% first get full timeseries values
[PSTHCorrs{1}, CorrTags] = ReplayCrossCorr(OpenLoopPSTH,reps_per_condition); % pair wise correlation of the single rep PSTHs
[PSTHResiduals{1}, ResidualTags] = ReplayResiduals(OpenLoopPSTH,reps_per_condition); % pair wise correlation of the single rep PSTHs

% Parse continuous traces into trials
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

% Split the long trace into odor-specific stretches 
% only keep points from this trial's odorstart to next trial's odor start
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

    [PSTHCorrs{1+whichodor}] = ReplayCrossCorr(OpenLoopPSTH(:,MyIdx,:),reps_per_condition);
    [PSTHResiduals{1+whichodor}] = ReplayResiduals(OpenLoopPSTH(:,MyIdx,:),reps_per_condition);
end


%% getting means and errors across pairs of either residuals or corrs
% correlations
U = unique(CorrTags); % various types of comparisons - Cl-OL, OL-OL, OL-PR etc
for i = 1:4
    for x = 1:length(U)
        CorrsMedian{i}(:,x) = median(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')';
        CorrsMean{i}(:,x)   = mean(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')';
        CorrsSTD{i}(:,x)    = std(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')';
        %ts = tinv([0.025  0.975],length(x)-1);      % T-Score
        nsamps              = numel(find(CorrTags==U(x)));
        ts                  = tinv([0.025  0.975],(nsamps-1));
        CorrsCI95{i}(:,x)   = ts(2)*std(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')'/sqrt(nsamps);
    end
end

% residuals
U = unique(ResidualTags); % various types of comparisons - Cl-OL, OL-OL, OL-PR etc
for i = 1:4
    for x = 1:length(U)
        ResidualsMedian{i}(:,x)     = median(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        ResidualsMean{i}(:,x)       = mean(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        ResidualsSTD{i}(:,x)        = std(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        nsamps                      = numel(find(ResidualTags==U(x)));
        ts                          = tinv([0.025  0.975],(nsamps-1));
        ResidualsCI95{i}(:,x)       = ts(2)*std(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')'/sqrt(nsamps);
    end
end
%% Plotting the outcomes
% Reorder Units by decreasing OL-OL correlation
[~,SortedbyCorr] = sort(CorrsMedian{1}(:,3),'descend');
SortedbyCorr = circshift(SortedbyCorr,-1); % first value was NaN'

whichtype = 1;
MedianCorrs = CorrsMedian{whichtype};
STDCorrs = CorrsSTD{whichtype};
MedianResiduals = ResidualsMean{whichtype};
STDResiduals = ResidualsSTD{whichtype};

figure;
subplot(4,1,1);
xpts = (1:4:4*N);
% OL-OL
bar(xpts, MedianCorrs(SortedbyCorr,3),'Facecolor',[1 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,3) + STDCorrs(SortedbyCorr,3))'; (MedianCorrs(SortedbyCorr,3) - STDCorrs(SortedbyCorr,3))'], ...
    'color', [1 0 0], 'Linewidth', 2);
%CL-OL
bar(xpts, MedianCorrs(SortedbyCorr,1),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,1) + STDCorrs(SortedbyCorr,1))'; (MedianCorrs(SortedbyCorr,1) - STDCorrs(SortedbyCorr,1))'], ...
    'color', [0 0 0],'Linewidth', 2);

subplot(4,1,2);
% PR-PR
xpts = (2:4:4*N);
bar(xpts, MedianCorrs(SortedbyCorr,5),'Facecolor',[0.4 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,5) + STDCorrs(SortedbyCorr,5))'; (MedianCorrs(SortedbyCorr,5) - STDCorrs(SortedbyCorr,5))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
% CL-PR
bar(xpts, MedianCorrs(SortedbyCorr,2),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,2) + STDCorrs(SortedbyCorr,2))'; (MedianCorrs(SortedbyCorr,2) - STDCorrs(SortedbyCorr,2))'], ...
    'color', [0 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyCorr)));
xtickangle(gca,90);

subplot(4,1,3);
xpts = (1:4:4*N);
% OL-OL
bar(xpts, MedianResiduals(SortedbyCorr,3),'Facecolor',[1 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,3) + STDResiduals(SortedbyCorr,3))'; (MedianResiduals(SortedbyCorr,3) - STDResiduals(SortedbyCorr,3))'], ...
    'color', [1 0 0], 'Linewidth', 2);
%CL-OL
bar(xpts, MedianResiduals(SortedbyCorr,1),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,1) + STDResiduals(SortedbyCorr,1))'; (MedianResiduals(SortedbyCorr,1) - STDResiduals(SortedbyCorr,1))'], ...
    'color', [0 0 0],'Linewidth', 2);

subplot(4,1,4);
% PR-PR
xpts = (2:4:4*N);
bar(xpts, MedianResiduals(SortedbyCorr,5),'Facecolor',[0.4 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,5) + STDResiduals(SortedbyCorr,5))'; (MedianResiduals(SortedbyCorr,5) - STDResiduals(SortedbyCorr,5))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
% CL-PR
bar(xpts, MedianResiduals(SortedbyCorr,2),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,2) + STDResiduals(SortedbyCorr,2))'; (MedianResiduals(SortedbyCorr,2) - STDResiduals(SortedbyCorr,2))'], ...
    'color', [0 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyCorr)));
xtickangle(gca,90);

%% scatter plot of self vs. across condition residuals
figure, 
subplot(1,2,1), hold on
subplot(1,2,2), hold on
for whichtype = 2:4
    %whichtype = 4;
    MedianCorrs = CorrsMedian{whichtype};
    STDCorrs = CorrsSTD{whichtype};
%     MedianResiduals = ResidualsMedian{whichtype};
%     STDResiduals = ResidualsSTD{whichtype};

    MedianResiduals = ResidualsMean{whichtype};
    STDResiduals = ResidualsCI95{whichtype};

    subplot(1,2,1)
    line([MedianResiduals(:,5) MedianResiduals(:,5)]', ...
        [(MedianResiduals(:,5) + STDResiduals(:,5)) (MedianResiduals(:,5) - STDResiduals(:,5))]', ...
        'Color', 'k');
    plot(MedianResiduals(:,5),MedianResiduals(:,2),'or');
%     line([MedianResiduals(:,5) MedianResiduals(:,5)]', ...
%         [(MedianResiduals(:,2) + STDResiduals(:,2)) (MedianResiduals(:,2) - STDResiduals(:,2))]', ...
%         'Color', 'r');
    
    subplot(1,2,2)
    line([MedianResiduals(:,3) MedianResiduals(:,3)]', ...
        [(MedianResiduals(:,3) + STDResiduals(:,3)) (MedianResiduals(:,3) - STDResiduals(:,3))]', ...
        'Color', 'k');
    plot(MedianResiduals(:,3),MedianResiduals(:,1),'or');
end

%%
figure;
N = length(MedianCorrs);
axes('ColorOrder',brewermap(N,'Spectral'),'NextPlot','replacechildren')

%
% plot within OL corrs on x axis
x1 = MedianCorrs(:,3) - [STDCorrs(:,3) -STDCorrs(:,3)];
y1 = [MedianCorrs(:,1) MedianCorrs(:,1)];
% plot OL-CL corr on Y axis
x2 = [MedianCorrs(:,3) MedianCorrs(:,3)];
y2 = MedianCorrs(:,1) - [STDCorrs(:,1) -STDCorrs(:,1)];
subplot(1,3,1)
line(x1',y1','Linewidth',1); %,'color',Plot_Colors('r'))
hold on
line(x2',y2','Linewidth',1); %,'color',Plot_Colors('r'))
line([-0.2 1], [-0.2 1], 'color', 'k', 'LineStyle',':');
set(gca,'XLim',[-0.2 1], 'YLim', [-0.2 1]);
axis square


% plot within PR corrs on x axis
x1 = MedianCorrs(:,5) - [STDCorrs(:,5) -STDCorrs(:,5)];
y1 = [MedianCorrs(:,2) MedianCorrs(:,2)];
% plot OL-CL corr on Y axis
x2 = [MedianCorrs(:,5) MedianCorrs(:,5)];
y2 = MedianCorrs(:,2) - [STDCorrs(:,2) -STDCorrs(:,2)];
subplot(1,3,2)
line(x1',y1','Linewidth',1); %,'color',Plot_Colors('t'))
hold on
line(x2',y2','Linewidth',1); %,'color',Plot_Colors('t'))
line([-0.2 1], [-0.2 1], 'color', 'k', 'LineStyle',':');
set(gca,'XLim',[-0.2 1], 'YLim', [-0.2 1]);
axis square

% plot within OL corrs on x axis
x1 = MedianCorrs(:,3) - [STDCorrs(:,3) -STDCorrs(:,3)];
y1 = [MedianCorrs(:,4) MedianCorrs(:,4)];
% plot OL-CL corr on Y axis
x2 = [MedianCorrs(:,3) MedianCorrs(:,3)];
y2 = MedianCorrs(:,4) - [STDCorrs(:,4) -STDCorrs(:,4)];
subplot(1,3,3)
line(x1',y1','Linewidth',1); %,'color',Plot_Colors('p'))
hold on
line(x2',y2','Linewidth',1); %,'color',Plot_Colors('p'))
line([-0.2 1], [-0.2 1], 'color', 'k', 'LineStyle',':');
set(gca,'XLim',[-0.2 1], 'YLim', [-0.2 1]);
axis square


%%
% Plot specific Units
%PlotUnits = [58 35 34 55 21];
PlotUnits = [34 55 27 32 30 51 6 26 57 52 41 46 12 44 28];
PlotUnits = [5 18 60 10 12];
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotephys', 1, 'UnitsPerFig', 5, 'whichunits', PlotUnits, 'PlotOpenLoop', 1);

