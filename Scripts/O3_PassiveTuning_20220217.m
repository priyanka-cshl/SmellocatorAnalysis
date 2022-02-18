%% for O3
SessionPath = 'O3/O3_20210927_r0_processed.mat';
ChosenUnits = [9 15 17 28 29]; %MyUnits = [8 35 28 55 39]; 

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

%% get trial aligned spikes for passive tuning
[TuningSpikes, TuningTTLs] = TrialAlignedSpikeTimesPassiveTuning(SingleUnits,TTLs,TuningTTLs);

%%
[SpikeCounts,Odors,Locations] = PlotPassiveTuning_v2(TuningSpikes, TuningTTLs); %, 'whichunits', [], 'rasters', 1, 'psth', 1);
% SpikeCounts - locations x [] x units x odors
% [noair air odor postodor]

%% get tuning curves for clopsed loop
[ClosedLoopCurve, XBins] = GetOdorTuningCurves(SessionPath,[],'binsize',50,'tuningbins',15);

%% Append CL and Passive tuning curves together
for i = 1:size(SingleUnits,2)
    TuningCurve(:,1:3,1,i) = squeeze(ClosedLoopCurve.ClosedLoopFull(:,1,i,:)); % Closed loop
    TuningCurve(:,1:3,2,i) = squeeze(SpikeCounts(:,3,i,2:4)); % Passive1
    TuningCurve(:,1:3,3,i) = squeeze(SpikeCounts(:,4,i,2:4)); % Passive1
end

%% Get Odor response values from TrialStart and from pertubation period

TrialInfo.TargetEntry = NaN*TrialInfo.Odor;
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

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
whichOdor = 1; % this was the perturbed odor
for i = 1:size(SingleUnits,2)
    % get spike counts aligned to odor ON
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 2);
    % count spikes in a 500 ms bin after odor ON
    for tz = 1:12
        bins2use = -BinOffset-500:-BinOffset;
        OdorResponse(tz,i) = mean(RawSpikeCounts(tz,bins2use)); 
    end
    % get spike counts aligned to perturbation start
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 6);
    % count spikes in a 500 ms bin after odor ON
    for tz = 1:12
        bins2use = -BinOffset:(-BinOffset+500);
        OdorResponse(tz,i) = mean(RawSpikeCounts(tz,bins2use)); 
    end
end

%% plot some examples
nrows = 4; ncols = 6;
for i = 1:size(SingleUnits,2)
    if mod(i,nrows*ncols) == 1
        figure;
        plotcount = 0;
    end
    plotcount = plotcount + 1;
    subplot(nrows,ncols,plotcount);
    hold on
    plot(mean(XBins,2),squeeze(TuningCurve(:,1,1,i)),'color',Plot_Colors('k'));
    plot(mean(XBins,2),squeeze(TuningCurve(:,1,3,i)),'color',Plot_Colors('t'));
    plot(mean(XBins,2),squeeze(TuningCurve(:,1,2,i)),'--','color',Plot_Colors('t'));
    title(num2str(i));
end