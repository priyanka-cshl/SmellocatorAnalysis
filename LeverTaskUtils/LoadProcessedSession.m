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

%% Concatenate traces and get one matrix with all behavior variables
% SampleRate = behavior sample rate;
[TracesOut, ColNames] = ConcatenateTraces2Mat(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);

%% sniff trace
%sniffs = TracesOut(:,3);
fband = [0.5 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
TracesOut(:,3) = filtfilt(b,a,TracesOut(:,3)); %apply the filter to x(t)

%% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior; % add to behavior TS to convert to OEPS

%% Process the TrialOn vector to add to it OdorOn periods in -ve
TrialTrace = TracesOut(:,find(strcmp(ColNames,'Trial')))';
% for every trial - patch in the odor ON period
for i = 1:size(TrialInfo.TrialID,2)
    x2 = TrialInfo.SessionIndices(i,1); % replay start
    x1 = x2 + TrialInfo.OdorStart(i,1)*SampleRate;
    if ~isnan(x1)
        TrialTrace(1,x1:x2) = -TrialInfo.Odor(i);
    end
end

%% If there are replay trials
if any(find(strcmp(TrialInfo.Perturbation,'OL-Replay')))
    OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);
    
    % for every trial - patch in the odor ON period
    replayTrials = OpenLoop.TemplateTraces.Trial{1};
    replayTrials(isnan(replayTrials)) = [];
    % only keep trace until last trial was off
    tracelength = find(diff(abs(replayTrials)==-1),1,'last');
    replayTrials(tracelength:end) = [];
    % remove initial samples - before trial ON
    replayTrials(1:find(replayTrials>0,1,'first')-1) = [];

    % patch in the trial sequence for the replay trials
    whichTrials = OpenLoop.ReplayTraces.TrialIDs{1};
    for i = 1:numel(whichTrials)
        x1 = TrialInfo.SessionIndices(whichTrials(i),1); % replay start
        x2 = x1 + numel(replayTrials) - 1;
        TrialTrace(1,x1:x2) = replayTrials*10;
    end
end

%% Split the trial vector into the three odors
OdorTrace = TrialTrace;
OdorTrace(abs(OdorTrace)>5) = OdorTrace(abs(OdorTrace)>5)/10;
for i = 1:3
    Motor = TracesOut(:,find(strcmp(ColNames,'Motor')))';
    Motor(abs(OdorTrace)~=i) = NaN;
    Motor = 125 - Motor;
    %Motor(isnan(Motor)) = 0;
    TracesOut(:,end+1) = Motor;
    ColNames{end+1} = ['Odor',num2str(i)];
end
% during air periods
Motor = TracesOut(:,find(strcmp(ColNames,'Motor')))';
Motor(abs(OdorTrace)~=0) = NaN;
Motor = 125 - Motor;
%Motor(isnan(Motor)) = 0;
TracesOut(:,end+1) = Motor;
ColNames{end+1} = ['Odor',num2str(4)];
