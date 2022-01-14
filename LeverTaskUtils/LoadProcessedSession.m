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