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
clear Traces
% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(ColNames,'Timestamps')))'; % in behavior timebase

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior; % add to behavior TS to convert to OEPS

%% Process the TrialOn vector to add to it OdorOn periods in -ve
TrialTrace = TracesOut(:,find(strcmp(ColNames,'Trial')))';
% get trial ON-OFF indices
Idx = [find(diff(TrialTrace>0)==1)'+1 find(diff(TrialTrace>0)==-1)'];
% get OdorStart Times w.r.t. Trial start (from the behavior file)
Idx(:,3) = Idx(:,1) + ceil(SampleRate*TrialInfo.OdorStart(whichTrials,2));
if ~isempty(TTLs)
    % OdorStart Times can also be extracted from OEPS TTLs
    Idx(:,4) = Idx(:,1) + ceil(SampleRate*TTLs.Trial(whichTrials,4));
    if any(abs(Idx(:,3)-Idx(:,4))>5)
        disp('mismatch in OdorOn Timestamps between behavior and OEPS files');
        keyboard;
    end
end
    temp = 0*TrialTrace;
    x1 = 1;
    for k = 1:size(Idx,1)
        % TargetZone vector
        x2 = Idx(k,2);
        temp(x1:x2,1) = AllTargets(TrialInfo.TargetZoneType(whichTrials(k)));
        x1 = x2 + 1;
        
        % OdorON + TrialON vector
        TrialTrace(Idx(k,3):Idx(k,1)-1,1) = -TrialTrace(Idx(k,1));
    end



%% Split the trial vector into the three odors
for i = 1:3
    Trial = TracesOut(:,find(strcmp(ColNames,'Trial')))';
    Trial(Trial~=i) = 0;
    Trial(Trial>0) = 1;
    Motor = TracesOut(:,find(strcmp(ColNames,'Motor')))';
    Motor(Trial==0) = NaN;
    Motor = 125 - Motor;
    Motor(Trial==0) = 0;
    TracesOut(:,end+1) = Motor;
    ColNames{end+1} = ['Odor',num2str(i)];
end