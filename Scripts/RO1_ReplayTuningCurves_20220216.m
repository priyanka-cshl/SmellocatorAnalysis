%%TuningCurve comparisons for Replay
%% for O3
SessionPath = 'O3/O3_20211005_r0_processed.mat';
ChosenUnits = [8 21 28 55 39]; %MyUnits = [8 35 28 55 39]; 

%% for PCX4
% SessionPath = 'PCX4/PCX4_20210721_r0_processed.mat';
% ChosenUnits = [2 38 45 55 64]; % may be 38

%% get tuning curves
% set binsize in SmellocatorTuning.m
% use fine bins (10) or coarse bins (24)
% fine bins give better tuning curves, but noisy residual comparisons

[TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves(SessionPath);

%% get cumulative distributions of the pairwise tuning curve correlations
XVar = (0:0.0001:0.35)';
AllResiduals = [PairedResiduals(:,1:3) ControlResiduals(:,1)];
[ResidualDist] = CumDist(AllResiduals,XVar);
cutoff = XVar(find(ResidualDist(:,4)>=0.95,1,'first'));

figure
hold on
axis square
area(XVar,ResidualDist(:,4),'LineStyle','none','FaceColor',[0.8 0.8 0.8]);
plot(XVar,ResidualDist(:,1),'Color',Plot_Colors('r'),'Linewidth',2);
plot(XVar,ResidualDist(:,2),'Color',Plot_Colors('t'),'Linewidth',2);
plot(XVar,ResidualDist(:,3),':','Color',Plot_Colors('t'),'Linewidth',2);
line([0 0.35],[0.95 0.95],'Linestyle',':','Color','k');
line(cutoff*[1 1], [0 1],'Linestyle',':','Color','k');

% Modulated units @95% cutoff
for i = 1:4
    counts(1,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff),4)));
end

% Modulated units @97% cutoff
cutoff97 = XVar(find(ResidualDist(:,4)>=0.97,1,'first'));
line([0 0.35],[0.97 0.97],'Linestyle',':','Color','b');
line(cutoff97*[1 1], [0 1],'Linestyle',':','Color','b');
for i = 1:4
    counts(2,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff97),4)));
end

% Modulated units @99% cutoff
cutoff99 = XVar(find(ResidualDist(:,4)>=0.99,1,'first'));
line([0 0.35],[0.99 0.99],'Linestyle',':','Color','r');
line(cutoff99*[1 1], [0 1],'Linestyle',':','Color','r');
for i = 1:4
    counts(3,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff99),4)));
end

title(mat2str([counts(1,:) counts(2,:) counts(3,:)]))

%% Plot the tuning curves
for unit = 1:numel(ChosenUnits)  %max(PairedResiduals(:,4))
    if mod(unit,5) == 1
        figure;
        i = 0;
    end
    i = i + 1;
    for j = 1:3
        subplot(5,3,j+3*(i-1));
        hold on;
        plot(mean(XBins,2)'-125,TuningCurve.ClosedLoopFull(:,1,ChosenUnits(unit),j),...
            'color',Plot_Colors('k'),'Linewidth',2);
        plot(mean(XBins,2)'-125,TuningCurve.OpenLoop(:,1,ChosenUnits(unit),j),...
            'color',Plot_Colors('r'),'Linewidth',2);
        plot(mean(XBins,2)'-125,TuningCurve.Passive(:,1,ChosenUnits(unit),j),...
            'color',Plot_Colors('t'),'Linewidth',2);
        
        foo = intersect(find(PairedResiduals(:,4)==ChosenUnits(unit)),find(PairedResiduals(:,5)==j));
        if any(PairedResiduals(foo,1:3)>cutoff)
            title(mat2str([ChosenUnits(unit) PairedResiduals(foo,1:3)],2),'Color','r');
        else
            title(mat2str([ChosenUnits(unit) PairedResiduals(foo,1:3)],2));
        end
    end
end

%% Plot the selected units
PlotReplayResponses(SessionPath,ChosenUnits);