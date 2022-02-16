% O3 - TuningCurve comparisons for Replay

%% get tuning curves
SessionPath = 'O3/O3_20211005_r0_processed.mat';
SessionPath = 'PCX4/PCX4_20210721_r0_processed.mat';
[TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves(SessionPath);

%% 
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

%% Modulated units/response
% @95%
for i = 1:4
    counts(1,i) = numel(unique(AllResiduals(find(AllResiduals(:,i)>cutoff),4)));
end

ChosenUnitsO3 = [8 21 28 55 39]; %MyUnits = [8 35 28 55 39]; 
ChosenUnitsPCX4 = [2 45 55 64]; % may be 38

%% 
for unit = 1:max(PairedResiduals(:,4))
    if mod(unit,5) == 1
        figure;
        i = 0;
    end
    i = i + 1;
    for j = 1:3
        subplot(5,3,j+3*(i-1));
        hold on;
        plot(mean(XBins,2)'-125,TuningCurve.ClosedLoopFull(:,1,unit,j),...
            'color',Plot_Colors('k'),'Linewidth',2);
        plot(mean(XBins,2)'-125,TuningCurve.OpenLoop(:,1,unit,j),...
            'color',Plot_Colors('r'),'Linewidth',2);
        plot(mean(XBins,2)'-125,TuningCurve.Passive(:,1,unit,j),...
            'color',Plot_Colors('t'),'Linewidth',2);
        
        foo = intersect(find(PairedResiduals(:,4)==unit),find(PairedResiduals(:,5)==j));
        if any(PairedResiduals(foo,1:3)>cutoff)
            title(mat2str([unit PairedResiduals(foo,1:3)],2),'Color','r');
        else
            title(mat2str([unit PairedResiduals(foo,1:3)],2));
        end
    end
end

%%
PlotReplayResponses(SessionPath);