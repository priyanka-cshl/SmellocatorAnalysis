%load '/mnt/data/DataMatrices/TuningCurves_O3_20211005.mat';
% Plotting Tuning Curves and noting residuals

%%
load('O3_20211005_TuningCurves_binsize10.mat');
Binsize = 10;
XBins = [];
XBins(:,1) = 20:Binsize:(230-Binsize);
XBins(:,2) = XBins(:,1) + Binsize;

%%
load('O3_20211005_TuningCurves_binsize24.mat');
Binsize = 24;
XBins = [];
XBins(:,1) = 24:Binsize:(224-Binsize);
XBins(:,2) = XBins(:,1) + Binsize;

%%
N = size(Curve_CL{1},4); % #units
% {} - different odors
% for each {}
% dimension 1 - bins
% dimension 2 - mean, median, sem, counts
% dimension 3 - actual, shuffled
% dimension 4 - units
Z = [];
for unit = 1:N
    for odor = 1:3
        X = [];
        a = max(Curve_CL{odor}(:,1,1,unit));
        X(:,1) = Curve_CL{odor}(:,1,1,unit)/a;
        X(:,2) = Curve_OL{odor}(:,1,1,unit)/a;
        X(:,3) = Curve_PR{odor}(:,1,1,unit)/a;
        X(:,4) = Curve_CL1{odor}(:,1,1,unit)/a;
        X(:,5) = Curve_CL2{odor}(:,1,1,unit)/a;
        
        if Binsize == 10
            X([1:4 21],:) = []; % delete unsmapled locations
        else
            X(1,:) = []; % delete unsmapled locations
        end
        
        if any(isnan(X))
            keyboard;
        end
        
        % calculate residuals
        Y(1) = mean((X(:,1) - X(:,2)).^2); % CL-OL
        Y(2) = mean((X(:,1) - X(:,3)).^2); % CL-PR
        Y(3) = mean((X(:,2) - X(:,3)).^2); % OL-PR
        Y(4) = mean((X(:,4) - X(:,5)).^2); % CL-CL
                
        Z = vertcat(Z, [Y' (1:4)' unit*ones(4,1) odor*ones(4,1)]);
    end
end

%% plotting cumulative distributions
figure;
% Z
% dimension 1 = residuals
% dimension 2 = paired condition for which the residual was calculated
% dimension 3 = unit identity
% dimension 4 = odor identity
XVar = (0:0.01:1)'; %logspace(-1,2)'; %
MyColors = [Plot_Colors('r'); Plot_Colors('t'); Plot_Colors('k')];
set(groot,'defaultAxesColorOrder',MyColors);
for odor = 1:4
    % binsize 10
    subplot(1,4,odor); hold on
    axis square
    if odor < 4
        AllResiduals = [];
        for conditions = 1:4
            AllResiduals(:,conditions) = Z(intersect(find(Z(:,2)==conditions),find(Z(:,4)==odor)),1);
        end
        [ResidualDist] = CumDist(AllResiduals,XVar);
    else
        AllResiduals = [];
        for conditions = 1:4
            AllResiduals(:,conditions) = Z(find(Z(:,2)==conditions),1);
        end
        [ResidualDist] = CumDist(AllResiduals,XVar);
    end
    area(XVar,ResidualDist(:,4),'LineStyle','none','FaceColor',[0.8 0.8 0.8]);
    plot(XVar,ResidualDist(:,1:2),'Linewidth',2);
    plot(XVar,ResidualDist(:,3),':','Color',Plot_Colors('t'),'Linewidth',2);
    line([0 1],[0.5 0.5],'Linestyle','--','Color','k');
end

%% get some numbers
cutoff = (find(ResidualDist(:,4)<=0.95,1,'last'));

MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat';

% get the processed data loaded and assemble open loop traces
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

ZSort = sortrows(Z,1,'descend');
ModUnits = unique(ZSort(1:42,3)); % top 10 residuals
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, 'plotfigures',1,'plotephys',1,'whichunits', ModUnits);

%% Lets get their Tuning Curves as well
figure;
for n = 1:numel(ModUnits) %size(Curve_CL{1},4)
    i = ModUnits(n);
    if mod(n,5) == 1
        figure;
        plotcount = 0;
    end
    for odor = 1:3
        plotcount = plotcount + 1;
        subplot(5,3,plotcount); hold on
        %line([0 0],[0 20],'LineStyle',':','Color','k');
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,1,i)',Curve_CL{odor}(:,3,1,i)',Plot_Colors('k'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,2,i)',Curve_CL{odor}(:,3,2,i)',Plot_Colors('k'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,1,i)',Curve_OL{odor}(:,3,1,i)',Plot_Colors('r'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,2,i)',Curve_OL{odor}(:,3,2,i)',Plot_Colors('r'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,1,i)',Curve_PR{odor}(:,3,1,i)',Plot_Colors('t'),{},0.5);
        %MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,2,i)',Curve_PR{odor}(:,3,2,i)',Plot_Colors('t'));
        
        r = Z(intersect(find(Z(:,3)==n),find(Z(:,4)==odor)),1:2);
        title(mat2str([i r(1:3,1)'],2));
    end
end

%%
ZSort = sortrows(Z,1,'ascend');
ModUnits = unique(ZSort(1:14,3));
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, 'plotfigures',1,'plotephys',1,'whichunits', ModUnits);
