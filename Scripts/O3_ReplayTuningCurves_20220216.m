% O3 - TuningCurve comparisons for Replay

%% get tuning curves
SessionPath = 'O3/O3_20211005_r0_processed.mat';
[TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves(SessionPath);

%% 
AllResiduals = [PairedResiduals(:,1:3) ControlResiduals(:,1)];
[ResidualDist] = CumDist(AllResiduals,XVar);
figure
hold on
area(XVar,ResidualDist(:,4),'LineStyle','none','FaceColor',[0.8 0.8 0.8]);
plot(XVar,ResidualDist(:,1:2),'Linewidth',2);
plot(XVar,ResidualDist(:,3),':','Color',Plot_Colors('t'),'Linewidth',2);
line([0 0.4],[0.95 0.95],'Linestyle','--','Color','k');

%% 
figure;
for j = 1:3 
    for i = 1:5
        subplot(5,3,j+3*(i-1));
        hold on;
        plot(mean(XBins,2)'-125,TuningCurve.ClosedLoopFull(:,1,i,j),...
            'color',Plot_Colors('k'),'Linewidth',2);
        plot(mean(XBins,2)'-125,TuningCurve.OpenLoop(:,1,i,j),...
            'color',Plot_Colors('r'),'Linewidth',2);
        plot(mean(XBins,2)'-125,TuningCurve.Passive(:,1,i,j),...
            'color',Plot_Colors('t'),'Linewidth',2);
    end
end