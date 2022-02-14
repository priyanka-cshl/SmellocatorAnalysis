% Plotting Tuning Curves and noting residuals
load('O3_20211005_TuningCurves_binsize10.mat');
Binsize = 10;
XBins = [];
XBins(:,1) = 20:Binsize:(230-Binsize);
XBins(:,2) = XBins(:,1) + Binsize;
%  
MyUnits = [58 35 34 55 21];

%%
figure;
for n = 1:numel(MyUnits); size(Curve_CL{1},4)
    i = MyUnits(n);
    if mod(n,5) == 1
        %figure;
        plotcount = 0;
    end
    for odor = 1:3
        plotcount = plotcount + 1;
        subplot(5,3,plotcount); hold on
        line([0 0],[0 20],'LineStyle',':','Color','k');
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,1,i)',Curve_CL{odor}(:,3,1,i)',Plot_Colors('k'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,2,i)',Curve_CL{odor}(:,3,2,i)',Plot_Colors('k'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,1,i)',Curve_OL{odor}(:,3,1,i)',Plot_Colors('r'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,2,i)',Curve_OL{odor}(:,3,2,i)',Plot_Colors('r'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,1,i)',Curve_PR{odor}(:,3,1,i)',Plot_Colors('t'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,2,i)',Curve_PR{odor}(:,3,2,i)',Plot_Colors('t'));
    end
end

load('O3_20211005_TuningCurves_binsize24.mat');
Binsize = 24;
XBins = [];
XBins(:,1) = 24:Binsize:(224-Binsize);
XBins(:,2) = XBins(:,1) + Binsize;

figure;
for n = 1:numel(MyUnits); size(Curve_CL{1},4)
    i = MyUnits(n);
    if mod(n,5) == 1
        %figure;
        plotcount = 0;
    end
    for odor = 1:3
        plotcount = plotcount + 1;
        subplot(5,3,plotcount); hold on
        line([0 0],[0 20],'LineStyle',':','Color','k');
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,1,i)',Curve_CL{odor}(:,3,1,i)',Plot_Colors('k'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,2,i)',Curve_CL{odor}(:,3,2,i)',Plot_Colors('k'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,1,i)',Curve_OL{odor}(:,3,1,i)',Plot_Colors('r'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,2,i)',Curve_OL{odor}(:,3,2,i)',Plot_Colors('r'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,1,i)',Curve_PR{odor}(:,3,1,i)',Plot_Colors('t'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,2,i)',Curve_PR{odor}(:,3,2,i)',Plot_Colors('t'));
    end
end