MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up
%MySession = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/O3/O3_20211005_r0_processed.mat';
%MySession = '/mnt/data/Processed/Behavior/PCX4/PCX4_20210721_r0_processed.mat'; % session path - leave empty to get browser pop up
LoadProcessedSession; % loads relevant variables

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

N = size(SingleUnits,2);
% sort units by tetrode - to match session viewer
clear foo
for i = 1:N
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id];
end
[~,SortedByTetrodes] = sort(foo(:,1));
UnitOrder = SortedByTetrodes;

[MyTraces,timestamps,PSTH,Raster] = ...
        ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'PSTHsmooth', 100, ...
        'plotfigures', 0, 'whichunits', UnitOrder);
    
% TraceNames = {'Lever' 'Motor' 'Sniffs' 'Licks' 'Rewards' 'Trial' 'TargetZone'};

%% Correlation analysis
[C, g] = ReplayCrossCorr(PSTH,[1 10 5]);

U = unique(g);
for x = 1:length(U)
    MedianCorrs(:,x) = median(C(find(strcmp(g,U{x})),:))';
    STDCorrs(:,x)    = std(C(find(strcmp(g,U{x})),:))';
end

%%
figure;
N = length(MedianCorrs);
axes('ColorOrder',brewermap(N,'Spectral'),'NextPlot','replacechildren')

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
% Reorder Units by decreasing OL-OL correlation
[~,SortedbyCorr] = sort(MedianCorrs(:,3),'descend');
SortedbyCorr = circshift(SortedbyCorr,-1); % first value was NaN'

figure;
subplot(2,1,1);
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
%PR-OL
xpts = xpts + 1.5;
plot(xpts, MedianCorrs(SortedbyCorr,4), 'o', ...
    'MarkerSize',6,'MarkerEdgecolor',[0.4 0.6 0.6],'MarkerFacecolor',[0.4 0.6 0.6]);
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,4) + STDCorrs(SortedbyCorr,4))'; (MedianCorrs(SortedbyCorr,4) - STDCorrs(SortedbyCorr,4))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(SortedByTetrodes(SortedbyCorr)));
xtickangle(gca,90);

subplot(2,1,2);
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
xticklabels(num2str(SortedByTetrodes(SortedbyCorr)));
xtickangle(gca,90);

%%
% Plot specific Units
%PlotUnits = SortedByTetrodes(SortedbyCorr(51:55));
PlotUnits = [58 35 34 55 21];
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotephys', 1, 'UnitsPerFig', 6, 'whichunits', PlotUnits);

