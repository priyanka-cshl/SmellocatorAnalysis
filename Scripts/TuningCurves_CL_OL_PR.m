MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up
%MySession = '/mnt/data/Processed/Behavior/PCX4/PCX4_20210721_r0_processed.mat';

% get the data loaded
[TracesOut, ColNames, TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, PassiveIdx, ...
    OpenLoop] = ...
    LoadProcessedDataSession(MySession);
% LoadProcessedSession; % loads relevant variables
% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(ColNames,'Timestamps')))'; % in behavior timebase

%% Get all PSTHs
% spike data - convert all to FRs - in same sampling rate as the behavior
N = size(SingleUnits,2);

AllFRs = NaN*zeros(size(TracesOut,1),N);
for i = 1:N
    whichUnit = i; %SortedByTetrodes(i);
    thisUnitspikes = SingleUnits(whichUnit).spikes'; % in OEPS time base 
    thisUnitspikes = thisUnitspikes - TimestampAdjuster; % in behavior time base
    %thisUnitFR = MakePSTH(thisUnitspikes, 0, PSTHWindow, 'downsample', SampleRate);
    thisUnitFR = MakePSTH_v2(thisUnitspikes, 0, 'binsize', 1, 'downsample', SampleRate);
    FRlength = min(size(TracesOut,1), numel(thisUnitFR));
    AllFRs(1:FRlength,i) = thisUnitFR(1:FRlength);
end

%% Passive replays
PassiveTimestamps = (1:size(PassiveTracesOut,1))/SampleRate;
PassiveFRs = NaN*zeros(size(PassiveTracesOut,1),N);
PassiveReps = size(PassiveIdx,1);
for i = 1:N
    whichUnit = i; %SortedByTetrodes(i);
    thisUnitspikes = SingleUnits(whichUnit).spikes'; % in OEPS time base
    ConcatSpikes = [];
    for j = 1:PassiveReps
        whichTrial = ReplayTTLs.TrialID(end-PassiveReps+j);
        TrialStart_OEPS = TTLs.Trial(ReplayTTLs.TrialID(end-PassiveReps+j),1);
        TrialStart_PassiveTrace = PassiveTimestamps(PassiveIdx(j,1));
        x1 = TrialStart_OEPS - 1;
        x2 = TTLs.Trial(ReplayTTLs.TrialID(end-PassiveReps+j),2) + 0.5;
        thisReplaySpikes = thisUnitspikes(((thisUnitspikes>=x1)&(thisUnitspikes<=x2)));
        thisReplaySpikes = thisReplaySpikes - TrialStart_OEPS + TrialStart_PassiveTrace;
        ConcatSpikes = horzcat(ConcatSpikes,thisReplaySpikes);
    end
    thisUnitFR = MakePSTH_v2(ConcatSpikes, 0, 'binsize', 1, 'downsample', SampleRate);
    FRlength = min(size(PassiveTracesOut,1), numel(thisUnitFR));
    PassiveFRs(1:FRlength,i) = thisUnitFR(1:FRlength);
end

%% behavioral variables vs. FR
% Separate out openloop blocks
% start index of the replay trial
OpenLoopIndices = [TrialInfo.SessionIndices(OpenLoop.ReplayTraces.TrialIDs{1},1) - mode(TrialInfo.TimeIndices(:,1)) ...
    TrialInfo.SessionIndices(OpenLoop.ReplayTraces.TrialIDs{1},2) + mode(TrialInfo.TimeIndices(:,1))];

Indices2exclude = 1+0*TracesOut(:,1);
for i = 1:size(OpenLoopIndices,1)
    Indices2exclude(OpenLoopIndices(i,1):OpenLoopIndices(i,2)) = 0;
end

%% work with few units
MyUnits = [58 35 34 55 21];
figure;
for i = 1:numel(MyUnits)
    
    % 1. Lever
    %scatter(TracesOut(find(Indices2exclude),1)',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    % 1. Sniffing
    %scatter(filteredsniffs(find(Indices2exclude))',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    % 2. Odor
    for odor = 1:4
        subplot(4,5,i+5*(odor-1)); hold on
        whichcol = 8+odor;
        scatter(TracesOut(find(Indices2exclude),whichcol)',AllFRs(find(Indices2exclude),MyUnits(i))',1);
        scatter(TracesOut(find(~Indices2exclude),whichcol)',AllFRs(find(~Indices2exclude),MyUnits(i))',0.5,'r');
        if odor < 4
            scatter(PassiveTracesOut(:,whichcol-1)',PassiveFRs(:,MyUnits(i))',0.5,Plot_Colors('t'));
        end
    end
end

%% Get Tuning Histograms - FR as a function of odor location
[Curve_CL, XBins] = SmellocatorTuning('Odor',TracesOut(find(Indices2exclude),9:12),AllFRs(find(Indices2exclude),MyUnits));
[Curve_OL, ~] = SmellocatorTuning('Odor',TracesOut(find(~Indices2exclude),9:12),AllFRs(find(~Indices2exclude),MyUnits));
[Curve_PR, ~] = SmellocatorTuning('Odor',PassiveTracesOut(:,8:11),PassiveFRs(:,MyUnits));

%%figure
%%
figure;
for i = 1:numel(MyUnits)
    for odor = 1:4
        subplot(4,5,i+5*(odor-1)); hold on
        MyShadedErrorBar(mean(XBins,2)',Curve_CL(:,1,odor,i)',Curve_CL(:,3,odor,i)',Plot_Colors('k'));
        MyShadedErrorBar(mean(XBins,2)',Curve_OL(:,1,odor,i)',Curve_OL(:,3,odor,i)',Plot_Colors('r'));
        if odor < 4
            MyShadedErrorBar(mean(XBins,2)',Curve_PR(:,1,odor,i)',Curve_PR(:,3,odor,i)',Plot_Colors('t'));
        end
    end
end

%% Get Tuning Histograms - odor location as a function of FR
[Curve_CL, XBins] = SmellocatorTuning('FR',AllFRs(find(Indices2exclude),MyUnits),TracesOut(find(Indices2exclude),9:12));
[Curve_OL, ~] = SmellocatorTuning('FR',AllFRs(find(~Indices2exclude),MyUnits),TracesOut(find(~Indices2exclude),9:12));
[Curve_PR, ~] = SmellocatorTuning('FR',PassiveFRs(:,MyUnits),PassiveTracesOut(:,8:11));

%%figure
%%
figure;
for i = 1:numel(MyUnits)
    for odor = 1:4
        subplot(4,5,i+5*(odor-1)); hold on
        MyShadedErrorBar(mean(XBins,2)',Curve_CL(:,1,odor,i)',Curve_CL(:,3,odor,i)',Plot_Colors('k'));
        MyShadedErrorBar(mean(XBins,2)',Curve_OL(:,1,odor,i)',Curve_OL(:,3,odor,i)',Plot_Colors('r'));
        if odor < 4
            MyShadedErrorBar(mean(XBins,2)',Curve_PR(:,1,odor,i)',Curve_PR(:,3,odor,i)',Plot_Colors('t'));
        end
    end
end


%%
figure;
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotephys', 1, 'UnitsPerFig', 5, 'whichunits', MyUnits);

    



