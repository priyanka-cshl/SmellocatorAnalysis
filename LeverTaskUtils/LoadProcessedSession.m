if ~exist('MySession','var')
    MySession = [];
end

%% Select session (if user didn't specify the path and load the data
if isempty(MySession)
    [Paths] = WhichComputer();
    [WhichSession, SessionPath] = uigetfile(...
        fullfile(Paths.ProcessedSessions,'O3/O3_20210922_r0_processed.mat'),...
        'Select Behavior or Recording Session');
    MySession = fullfile(SessionPath,WhichSession);
end

% Load the relevant variables
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

%SessionLength = 10*ceil(TTLs.Trial(end,2)/10);
%NumUnits = size(SingleUnits,2);
if find(strcmp(TrialInfo.Perturbation,'OL-Replay'));
OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);
end

%% Concatenate traces and get one matrix with all behavior variables
% SampleRate = behavior sample rate;
[TracesOut, ColNames] = ConcatenateTraces2Mat(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);
% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(ColNames,'Timestamps')))'; % in behavior timebase

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior; % add to behavior TS to convert to OEPS

%% Process the TrialOn vector to add to it OdorOn periods in -ve
TrialTrace = TracesOut(:,find(strcmp(ColNames,'Trial')))';
% for every trial - patch in the odor ON period
for i = 1:size(TrialInfo.TrialID,2)
    x2 = TrialInfo.SessionIndices(i,1);
    x1 = x2 + TrialInfo.OdorStart(i,2)*SampleRate;
    if ~isnan(x1)
        TrialTrace(1,x1:x2) = -TrialInfo.Odor(i);
    end
end

%% Split the trial vector into the three odors
for i = 1:3
    Motor = TracesOut(:,find(strcmp(ColNames,'Motor')))';
    Motor(abs(TrialTrace)~=i) = NaN;
    Motor = 125 - Motor;
    %Motor(isnan(Motor)) = 0;
    TracesOut(:,end+1) = Motor;
    ColNames{end+1} = ['Odor',num2str(i)];
end