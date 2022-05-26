% --- Executes on button press in LoadSession.
%function LeverTaskLoadSession
% hObject    handle to LoadSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[WhichSession, SessionPath] = uigetfile(...
                                fullfile('O3/O3_20210922_r0_processed.mat'),...
                                'Select Behavior or Recording Session');
WhereSession = fullfile(SessionPath,WhichSession);

% Load the relevant variables
load(WhereSession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

if ~isempty(TTLs)
    SessionLength = num2str(10*ceil(TTLs.Trial(end,2)/10));
    NumUnits = num2str(size(SingleUnits,2));
else
    SessionLength = TrialInfo.SessionTimestamps(end,2);
    NumUnits = 'NaN';
end

if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    u = unique(TrialInfo.Perturbation(x));
    PerturbationList = u{1};
    for y = 2:size(u,1)
        PerturbationList = [PerturbationList,'; ',u{y}];
    end
else
    x = [];
    PerturbationList = '';
end

%% Concatenate traces
whichTrials = 1:length(TrialInfo.TrialID);
traceOverlap = SampleRate*startoffset;
%whichTraces = fieldnames(Traces);
whichTraces{1} = 'Lever'; whichTraces{2} = 'Motor'; whichTraces{3} = 'Sniffs';
whichTraces{4} = 'Trial'; whichTraces{5} = 'Rewards'; whichTraces{6} = 'Timestamps';

for j = 1:size(whichTraces,2)
    temp = cellfun(@(x) ...
        x(1:end-traceOverlap), Traces.(whichTraces{j})(whichTrials), ...
        'UniformOutput', false);
    TracesOut.(whichTraces{j}) = {[cell2mat(temp(:)); ...
        Traces.(whichTraces{j}){whichTrials(end)}(end-traceOverlap+1:end,1)]};
end

% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);

if ~isempty(TTLs)
    TrialStart_Ephys = TTLs.Trial(1,2);
    % factor to convert all behavior timestamps to match Ephys
    TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;
else
    TrialStart_Ephys = 0;
    TimestampAdjuster = 0;
end

TracesOut.Timestamps{1} = Timestamps + TimestampAdjuster;

% Get OpenLoop data to get the Trial vector
[OpenLoop] = ProcessReplayResponses(WhereSession);

%% flag replay periods 
% find perturbation periods
PerturbTS = [];
if any(strcmp(TrialInfo.Perturbation(x),'OL-Template'))
    % get start and stop TS of the template
    templateStart = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'first'));
    templateEnd   = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'last'));
    PerturbTS(1,1) = TrialInfo.SessionTimestamps(templateStart,1) + TimestampAdjuster;
    PerturbTS(2,1) = TrialInfo.SessionTimestamps(templateEnd,2) + TimestampAdjuster;
    % get start and stop of the replays
    replays = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Replay')));
    PerturbTS = horzcat(PerturbTS, ...
        TrialInfo.SessionTimestamps(replays,1:2)' + TimestampAdjuster );
end

OpenLoopTrialVec = OpenLoop.Traces.Trial(:,1);
OpenLoopTrialVec(1:(find(OpenLoopTrialVec>0,1,'first')-1),:) = [];
OpenLoopTrialVec(OpenLoopTrialVec<0) = 0;

TracesOut.OpenLoop = TracesOut.Trial;
for i = 1:size(PerturbTS,2)
    idx1 = find(TracesOut.Timestamps{1}>=PerturbTS(1,i),1,'first') + 1;
    idx2 = find(TracesOut.Timestamps{1}<=PerturbTS(2,i),1,'last') - 1;
    idx3 = idx2 - idx1 + 1;
    
    if i > 1
        if (unique(TracesOut.Trial{1}(idx1:idx2)))
            %TracesOut.Trial{1}(idx1:idx2) = -1;
            TracesOut.OpenLoop{1}(idx1:idx2) = -2;
            TracesOut.Trial{1}(idx1:idx2) = OpenLoopTrialVec(1:idx3);
        else
            keyboard;
        end
    else
        TracesOut.OpenLoop{1}(idx1:idx2) = -1;
    end
end

temp = TracesOut.OpenLoop{1};
temp(temp>0) = 0;
TracesOut.OpenLoop{1} = abs(temp)';

%% Organize Open loop data separately

path = fileparts(fileparts(fileparts(WhereSession)));
if ~exist(fullfile(path,'Alex'),'dir')
    mkdir(fullfile(path,'Alex'));
end
[~,filename] = fileparts(WhereSession);
filename = [filename,'.mat'];
save(fullfile(path,'Alex',filename),'TracesOut','SingleUnits');

% clearvars -except TracesOut SingleUnits