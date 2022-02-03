MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up

% get the data loaded 
LoadProcessedSession; % loads relevant variables
% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(ColNames,'Timestamps')))'; % in behavior timebase

%% Get all PSTHs
% spike data - convert all to FRs - in same sampling rate as the behavior
N = size(SingleUnits,2);
% sort units by tetrode - to match session viewer
clear foo
for i = 1:N
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id];
end
[~,SortedByTetrodes] = sort(foo(:,1));

SessionLength = 10*ceil(Timestamps(end)/10); % last timestamp in the behavior file
PSTHWindow = [0 1000*SessionLength];
AllFRs = NaN*zeros(size(TracesOut,1),N);
for i = 1:N
    whichUnit = i; %SortedByTetrodes(i);
    thisUnitspikes = SingleUnits(whichUnit).spikes'; % in OEPS time base 
    thisUnitspikes = thisUnitspikes - TimestampAdjuster; % in behavior time base
    thisUnitFR = MakePSTH(thisUnitspikes, 0, PSTHWindow, 'downsample', SampleRate);
    AllFRs(:,i) = thisUnitFR(1:size(TracesOut,1));
end

%% Passive replays
LoadProcessedPassiveSession;
PassiveTimestamps = (1:size(PassiveTracesOut,1))/SampleRate;
PassiveFRs = NaN*zeros(size(PassiveTracesOut,1),N);
SessionLength = 10*ceil(PassiveTimestamps(end)/10); % last timestamp in the behavior file
PSTHWindow = [0 1000*SessionLength];
for i = 1:N
    whichUnit = i; %SortedByTetrodes(i);
    thisUnitspikes = SingleUnits(whichUnit).spikes'; % in OEPS time base
    ConcatSpikes = [];
    for j = 1:nReps
        whichTrial = ReplayTTLs.TrialID(end-nReps+j);
        TrialStart_OEPS = TTLs.Trial(ReplayTTLs.TrialID(end-nReps+j),1);
        TrialStart_PassiveTrace = PassiveTimestamps(StartStopIdx(j,1));
        x1 = TrialStart_OEPS - 1;
        x2 = TTLs.Trial(ReplayTTLs.TrialID(end-nReps+j),2) + 0.5;
        thisReplaySpikes = thisUnitspikes(((thisUnitspikes>=x1)&(thisUnitspikes<=x2)));
        thisReplaySpikes = thisReplaySpikes - TrialStart_OEPS + TrialStart_PassiveTrace;
        ConcatSpikes = horzcat(ConcatSpikes,thisReplaySpikes);
    end
    thisUnitFR = MakePSTH(ConcatSpikes, 0, PSTHWindow, 'downsample', SampleRate);
    PassiveFRs(:,i) = thisUnitFR(1:size(PassiveTracesOut,1));
end

%% behavioral variables vs. FR
% Separate out openloop blocks
% start index of the replay trial
OpenLoopIndices = [TrialInfo.SessionIndices(OpenLoop.ReplayTraces.TrialIDs{1},1) - startoffset*SampleRate ...
    TrialInfo.SessionIndices(OpenLoop.ReplayTraces.TrialIDs{1},2) + startoffset*SampleRate];

Indices2exclude = 1+0*Timestamps;
for i = 1:size(OpenLoopIndices,1)
    Indices2exclude(OpenLoopIndices(i,1):OpenLoopIndices(i,2)) = 0;
end

%% work with few units
MyUnits = [58 35 34 55 21];
figure;
for i = 1:numel(MyUnits)
    subplot(4,5,i); hold on
    % 1. Lever
    %scatter(TracesOut(find(Indices2exclude),1)',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    % 1. Sniffing
    %scatter(filteredsniffs(find(Indices2exclude))',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    % 2. Odor
    scatter(TracesOut(find(Indices2exclude),11)',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    scatter(TracesOut(find(~Indices2exclude),11)',AllFRs(find(~Indices2exclude),MyUnits(i))',0.5,'r');
    scatter(PassiveTracesOut(:,11)',PassiveFRs(:,MyUnits(i))',0.5,Plot_Colors('t'));
    subplot(4,5,i+5); hold on
    scatter(TracesOut(find(Indices2exclude),9)',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    scatter(TracesOut(find(~Indices2exclude),9)',AllFRs(find(~Indices2exclude),MyUnits(i))',0.5,'r');
    scatter(PassiveTracesOut(:,9)',PassiveFRs(:,MyUnits(i))',0.5,Plot_Colors('t'));
    subplot(4,5,i+10); hold on
    scatter(TracesOut(find(Indices2exclude),8)',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    scatter(TracesOut(find(~Indices2exclude),8)',AllFRs(find(~Indices2exclude),MyUnits(i))',0.5,'r');
    scatter(PassiveTracesOut(:,8)',PassiveFRs(:,MyUnits(i))',0.5,Plot_Colors('t'));
    subplot(4,5,i+15); hold on
    scatter(TracesOut(find(Indices2exclude),10)',AllFRs(find(Indices2exclude),MyUnits(i))',1);
    scatter(TracesOut(find(~Indices2exclude),10)',AllFRs(find(~Indices2exclude),MyUnits(i))',0.5,'r');
    scatter(PassiveTracesOut(:,10)',PassiveFRs(:,MyUnits(i))',0.5,Plot_Colors('t'));
end

%% Get Tuning Histograms
[Curve_CL, XBins] = SmellocatorTuning('Odor',TracesOut(find(Indices2exclude),8:11),AllFRs(find(Indices2exclude),MyUnits));
[Curve_OL, ~] = SmellocatorTuning('Odor',TracesOut(find(~Indices2exclude),8:11),AllFRs(find(~Indices2exclude),MyUnits));
[Curve_PR, ~] = SmellocatorTuning('Odor',PassiveTracesOut(:,8:11),PassiveFRs(:,MyUnits));

%%figure
%%
figure;
for i = 1:numel(MyUnits)
    for odor = 1:4
        subplot(4,5,i+5*(odor-1)); hold on
%         plot(mean(XBins,2),squeeze(Curve_CL(:,1,odor,i)),'color',Plot_Colors('b'));
%         plot(mean(XBins,2),squeeze(Curve_OL(:,1,odor,i)),'color',Plot_Colors('r'));
%         plot(mean(XBins,2),squeeze(Curve_PR(:,1,odor,i)),'color',Plot_Colors('t'));
        MyShadedErrorBar(mean(XBins,2)',Curve_CL(:,1,odor,i)',Curve_CL(:,3,odor,i)',Plot_Colors('k'));
        MyShadedErrorBar(mean(XBins,2)',Curve_OL(:,1,odor,i)',Curve_OL(:,3,odor,i)',Plot_Colors('r'));
        MyShadedErrorBar(mean(XBins,2)',Curve_PR(:,1,odor,i)',Curve_PR(:,3,odor,i)',Plot_Colors('t'));
    end
end

%%
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotephys', 1, 'UnitsPerFig', 5, 'whichunits', MyUnits);

    



