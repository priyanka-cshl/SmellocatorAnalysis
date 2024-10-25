function [] = PlottingSelectedReplayUnits(MyUnits, unitsperfig)

if nargin < 2
    unitsperfig = 5;
end

MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat';
MySession = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/O3/O3_20211005_r0_processed.mat';

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');
            
OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'PSTHsmooth', 100, ...
        'plotfigures', 0, 'plotephys', 1, 'UnitsPerFig', unitsperfig, ...
        'whichunits', MyUnits);

%% Also plot their Tuning curves
load ('/mnt/data/Processed/Behavior/O3/O3_20211005_TuningCurves_binsize10.mat');
Binsize = 10;
XBins = [];
XBins(:,1) = 20:Binsize:(230-Binsize);
XBins(:,2) = XBins(:,1) + Binsize;

figure;
for n = 1:numel(MyUnits) %size(Curve_CL{1},4)
    i = MyUnits(n);
    if mod(n,unitsperfig) == 1
        figure;
        whichrow = 0;
    end
    
    for odor = 1:3
        subplot(unitsperfig,3,(3*whichrow)+odor); 
        hold on
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL_10{odor}(:,1,1,i)',Curve_CL_10{odor}(:,3,1,i)',Plot_Colors('k'),{},0.5);
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL_10{odor}(:,1,1,i)',Curve_OL_10{odor}(:,3,1,i)',Plot_Colors('r'),{},0.5);
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR_10{odor}(:,1,1,i)',Curve_PR_10{odor}(:,3,1,i)',Plot_Colors('t'),{},0.5);
        
        r = Z10(intersect(find(Z10(:,3)==n),find(Z10(:,4)==odor)),1:2);
        
        myCurves = [Curve_CL_10{odor}(:,1,1,i),Curve_OL_10{odor}(:,1,1,i),Curve_PR_10{odor}(:,1,1,i)];
        foo = corrcoef(myCurves(5:20,:));
        
        title(mat2str([i r(1:3,1)' foo(2:3,1)' foo(3,2)],2));
    end
    whichrow = whichrow + 1;
end

%% Also plot their Tuning curves
load ('/mnt/data/Processed/Behavior/O3/O3_20211005_TuningCurves_binsize24.mat');
Binsize = 24;
XBins = [];
XBins(:,1) = 24:Binsize:(224-Binsize);
XBins(:,2) = XBins(:,1) + Binsize;

figure;
for n = 1:numel(MyUnits) %size(Curve_CL{1},4)
    i = MyUnits(n);
    if mod(n,unitsperfig) == 1
        figure;
        whichrow = 0;
    end
    
    for odor = 1:3
        subplot(unitsperfig,3,(3*whichrow)+odor); 
        hold on
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,1,i)',Curve_CL{odor}(:,3,1,i)',Plot_Colors('k'),{},0.5);
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,1,i)',Curve_OL{odor}(:,3,1,i)',Plot_Colors('r'),{},0.5);
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,1,i)',Curve_PR{odor}(:,3,1,i)',Plot_Colors('t'),{},0.5);
        
        r = Z24(intersect(find(Z24(:,3)==n),find(Z24(:,4)==odor)),1:2);
        
        myCurves = [Curve_CL{odor}(:,1,1,i),Curve_OL{odor}(:,1,1,i),Curve_PR{odor}(:,1,1,i)];
        foo = corrcoef(myCurves(2:end,:));
        
        title(mat2str([i r(1:3,1)' foo(2:3,1)' foo(3,2)],2));
    end
    whichrow = whichrow + 1;
end
%%
% O3 - 20211005
% UnModulated Units: [8 21 30 32 34 58] % 5 7 9 21 27 32 34 58
% [5 7 8 9 21 27 30 32 34 58]
% [8 9 27 35]
% Modulated Units:   [16 28 30 37 39 48 49 51 55]