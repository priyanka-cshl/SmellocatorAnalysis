%% for O3
SessionPath = 'O3/O3_20210927_r0_processed.mat';

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
            
[~, ~, TrialInfo] = LoadProcessedDataSession(MySession); % to get target zone entry time

%% Trial Aligned spikeTimes
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

%%
whichOdor = 1; 
for tz = 1:12
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
                    find(TrialInfo.Odor==whichOdor));
    whichTrials = intersect(whichTrials,...
                    find(TrialInfo.TargetZoneType==tz));
    OdorONPeriod(tz,1) = median(TrialInfo.OdorStart(whichTrials,1));
    foo = cellfun(@(x) median(x(491:500)), Traces.Motor(whichTrials),'UniformOutput', false);
    LocationTrialStart(tz,1) = round(median(cell2mat(foo)));
end

%% Plot
figure;
AlignTo = 6; 
switch AlignTo
    case {1,2,6}
        myXlim = [-1.2 6];
    case {3,4}
        myXlim = [-5.2 1];
    case 5
        myXlim = [-2.2 5];
end

whichUnit = 15; 
% whichTZ = 1;
trialsdone = 0;
subplot(3,2,[1 3]); hold on
for whichTZ = 2:2:12
trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
[nTrials, FRs, BinOffset] = ...
    UnitPlotter(whichUnit, whichOdor, whichTZ, AlignedSpikes, Events, TrialInfo, AlignTo, trialsdone, trialboxcolor);
trialsdone = trialsdone + nTrials;
end

% halt flip trials
for whichTZ = 1:12
trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
[nTrials, FRs, BinOffset] = ...
    UnitPlotter(whichUnit, whichOdor, -whichTZ, AlignedSpikes, Events, TrialInfo, AlignTo, trialsdone, trialboxcolor);
trialsdone = trialsdone + nTrials;
end
set(gca,'XLim',myXlim,'TickDir','out');


[FRs, BinOffset, whichTZ] = ... 
PlotHaltFlips(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo);
set(gca,'XLim',myXlim);

subplot(3,2,[5]); hold on
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);

whichUnit = 1;
subplot(3,2,[2 4]); hold on
[SRs, BinOffset, whichTZ] = ... 
PlotHaltFlips(whichUnit, whichOdor, AlignedSniffs, Events, TrialInfo, AlignTo);
set(gca,'XLim',myXlim);

subplot(3,2,[6]); hold on
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);
