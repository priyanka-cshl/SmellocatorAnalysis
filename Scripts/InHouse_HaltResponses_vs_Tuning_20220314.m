%% for O3
SessionPath = 'PCX4/PCX4_20210713_r0_processed.mat';

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

ChosenUnits = 1:size(SingleUnits,2); 

%% get the closed loop tuning curve
[TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves(SessionPath, ChosenUnits, 'tuningbins', 15); %, 'binsize', 4);

%% Get Odor response values from TrialStart and from pertubation period
TrialInfo.TargetEntry = NaN*TrialInfo.Odor;
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,1),TrialInfo,MySession);

whichOdor = 1;
%%
for tz = 1:12
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
                    find(TrialInfo.Odor==whichOdor));
    whichTrials = intersect(whichTrials,...
                    find(TrialInfo.TargetZoneType==tz));
    OdorONPeriod(tz,1) = median(TrialInfo.OdorStart(whichTrials,1));
    foo = cellfun(@(x) median(x(491:500)), Traces.Motor(whichTrials),'UniformOutput', false);
    LocationTrialStart(tz,1) = round(median(cell2mat(foo)));
end
% convert to Bins
OdorONPeriod = ceil(OdorONPeriod*SampleRate); 

% Align to:
% 1 - Trial ON, 2 - Odor ON, 3 - Trial OFF, 4 - Reward, 5 - 1st TZ entry, 6 - Perturbation Start
mywin = 1000; % in m
stepsize = 1000/SampleRate;
for i = 1:size(SingleUnits,2)
    % get spike counts aligned to odor ON
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 2);
    % count spikes in a 500 ms bin after odor ON
    bins = -BinOffset + [stepsize:stepsize:mywin];
    bins = bins/stepsize;
    foo = round(numel(bins)/2);
    for tz = 1:12
        % get odor response during 'windowsize' after odor ON
        AreaUnderCurve.Odor(tz,1,i) = sum(AlignedFRs(tz,bins));
        AreaUnderCurve.Odor(tz,2,i) = sum(AlignedFRs(tz,bins(1:foo)));
        AreaUnderCurve.Odor(tz,3,i) = sum(AlignedFRs(tz,bins((foo+1):end)));
    end
    % get spike counts aligned to perturbation start
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 6);
    % count spikes in a 500 ms bin after odor ON
    bins = -BinOffset + [stepsize:stepsize:mywin];
    bins = bins/stepsize;
    foo = round(numel(bins)/2);
    if size(AlignedPerturbationFRs,2)<bins(end)
        AlignedPerturbationFRs = horzcat(AlignedPerturbationFRs,...
            zeros(size(AlignedPerturbationFRs,1),(bins(end)- size(AlignedPerturbationFRs,2))));
    end
    for tz = 1:12
        AreaUnderCurve.Halt(tz,1,i) = sum(AlignedPerturbationFRs(tz,bins));
        AreaUnderCurve.Halt(tz,2,i) = sum(AlignedPerturbationFRs(tz,bins(1:foo)));
        AreaUnderCurve.Halt(tz,3,i) = sum(AlignedPerturbationFRs(tz,bins((foo+1):end)));
    end
end

AreaUnderCurve.Odor(:,1,:) = AreaUnderCurve.Odor(:,1,:)/(mywin/stepsize);
AreaUnderCurve.Halt(:,1,:) = AreaUnderCurve.Halt(:,1,:)/(mywin/stepsize);
AreaUnderCurve.Odor(:,2:3,:) = AreaUnderCurve.Odor(:,2:3,:)/(mywin/stepsize/2);
AreaUnderCurve.Halt(:,2:3,:) = AreaUnderCurve.Halt(:,2:3,:)/(mywin/stepsize/2);

%% Plot the effect
figure;
subplot(1,2,1);
hold on
axis square
subplot(1,2,2); 
hold on
axis square
whichcol = 1;
for i = 1:size(SingleUnits,2)
    % mean perturbation response vs. mean mirror location response
    subplot(1,2,1);
    plot(mean(AreaUnderCurve.Odor(11:12,whichcol,i)),mean(AreaUnderCurve.Halt([1 5 9],whichcol,i)),...
        'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
    
    line([0 40],[0 40],'Color','k')
    subplot(1,2,2);
    plot(max(AreaUnderCurve.Odor(:,whichcol,i)),mean(AreaUnderCurve.Halt([1 5 9],whichcol,i)),...
        'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
    line([0 40],[0 40],'Color','k')
end

%% halt response
% expected vs actual
HaltLocation = TrialInfo.Perturbation{find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)),1,'first'),2}(3);
WhichBin = find(mean(XBins,2)==HaltLocation);
Expected = squeeze(TuningCurve.ClosedLoopFull(WhichBin,:,:,whichOdor))';
Actual = squeeze(mean(AreaUnderCurve.Halt([1 5 9],1,:),1));
for i = 1:size(SingleUnits,2)
    Scored(i,1) = (Actual(i,1) - Expected(i,1))/Expected(i,2);
end
figure; hold on; axis square
plot(Expected(:,1),Actual(:,1), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
set(gca,'XLim',[0 40],'YLim',[0 40])
line([0 40],[0 40],'Color','k')

% ModUnits    = [10 15 17 22 28 29 31 32 52 63 66 71 73];
% UnModUnits  = [12 18 56];
% plot(Expected(ModUnits,1),Actual(ModUnits,1), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
% plot(Expected(UnModUnits,1),Actual(UnModUnits,1), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b');

%% Tuning curves vs. Odor On reponse
figure;
for i = 1:size(SingleUnits,2)
    subplot(7,11,i); hold on
    % CL
    plot(mean(XBins,2),squeeze(TuningCurve.ClosedLoopFull(:,1,i,whichOdor)),'k','Linewidth',2);
    % Passive
    %plot(Locations,squeeze(SpikeCounts(:,4,i,whichOdor+1)),'b','Linewidth',2);
    % Odor ON
    %plot(LocationTrialStart,AreaUnderCurve.Odor(:,1,i),'r','Linewidth',2);
    %plot(LocationTrialStart,AreaUnderCurve.Odor(:,2,i),'g','Linewidth',1);
    
    % Halt response
    plot(30,mean(AreaUnderCurve.Halt([1 5 9],1,i)),'.r');
end