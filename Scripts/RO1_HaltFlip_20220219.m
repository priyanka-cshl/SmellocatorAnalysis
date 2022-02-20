%% for O3
SessionPath = 'O3/O3_20210927_r0_processed.mat';
ChosenUnits = [15 18 28 63 56]; %MyUnits = [8 35 28 55 39]; 

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

%% Get Odor response values from TrialStart and from pertubation period

TrialInfo.TargetEntry = NaN*TrialInfo.Odor;
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);


whichUnit = 15;
whichOdor = 1;
AlignTo = 6; % perturbation start
MyColors1 = brewermap(15,'*PuBu');
figure;
for i = 1:numel(ChosenUnits)
    whichUnit = ChosenUnits(i);
    subplot(3,5,i); hold on
    [myFRs, BinOffset, whichTZ] = PlotHaltFlips(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo);
    set(gca, 'XLim', [-1.2 6]);
    title(['unit# ',num2str(whichUnit)]);
    
    subplot(3,5,i+5); hold on
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
    
    subplot(3,5,i+10); hold on
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
    % get FRs aligned to odor start as well - 
    [~, FRs, BinOffset] = PlotFullSession(-whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, 2);
    for t = 1:size(FRs,1)
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(t,:),'Color',MyColors1(t,:),'Linewidth',1);
    end
    set(gca, 'XLim', [-1.2 6]);
    
    % add odor start response to previous plot as well
    subplot(3,5,i+5); hold on
    plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(whichTZ,:),'color',Plot_Colors('r'),'Linewidth',1);
    set(gca, 'XLim', [-1.2 6]);
    
    
end

for tz = 1:12
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
                    find(TrialInfo.Odor==1));
    whichTrials = intersect(whichTrials,...
                    find(TrialInfo.TargetZoneType==tz));
    OdorONPeriod(tz,1) = median(TrialInfo.OdorStart(whichTrials,1));
    foo = cellfun(@(x) median(x(491:500)), Traces.Motor(whichTrials),'UniformOutput', false);
    LocationTrialStart(tz,1) = round(median(cell2mat(foo)));
end
% convert to Bins
OdorONPeriod = ceil(OdorONPeriod*SampleRate); 

% Align to:
% 1 - Trial ON, 2 - Odor ON, 3 - Trial OFF, 4 - Reward, 5 - 1st TZ entry, 
% 6 - Perturbation Start
%%
whichOdor = 1; % this was the perturbed odor
mywin = 1000; % in m
stepsize = 1000/SampleRate;
for i = 1:size(SingleUnits,2)
    % get spike counts aligned to odor ON
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 2);
    % count spikes in a 500 ms bin after odor ON
    bins = -BinOffset + [stepsize:stepsize:mywin];
    bins = bins/stepsize;
    for tz = 1:12
        % get odor response during 'windowsize' after odor ON
        AreaUnderCurve(tz,1,i) = sum(AlignedFRs(tz,bins));
    end
    % get spike counts aligned to perturbation start
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 6);
    % count spikes in a 500 ms bin after odor ON
    bins = -BinOffset + [stepsize:stepsize:mywin];
    bins = bins/stepsize;
    for tz = 1:12
        AreaUnderCurve(tz,2,i) = sum(AlignedPerturbationFRs(tz,bins));
    end
end
AreaUnderCurve = AreaUnderCurve/(mywin/stepsize);

%% Plot the effect
figure;
subplot(1,2,1); 
hold on
axis square
subplot(1,2,2); 
hold on
axis square
for i = 1:size(SingleUnits,2)
    % mean perturbation response vs. mean mirror location response
    subplot(1,2,1);
    plot(mean(AreaUnderCurve(11:12,1,i)),mean(AreaUnderCurve([1 5 9],2,i)),'ok');
    subplot(1,2,2);
    for tz = 1:12
        plot(mean(AreaUnderCurve(tz,1,i)),mean(AreaUnderCurve(1,2,i)),'ok');
        plot(mean(AreaUnderCurve(tz,1,i)),mean(AreaUnderCurve(5,2,i)),'or');
        plot(mean(AreaUnderCurve(tz,1,i)),mean(AreaUnderCurve(9,2,i)),'ob');
    end
end


%% plot some examples