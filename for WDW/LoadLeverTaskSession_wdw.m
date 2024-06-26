
% load the processed behavior/ephys files (from PreprocessSmellocatorData.m)
% [WhichSession, SessionPath] = uigetfile(...
%                                 fullfile('Q4/Q4_20221109_r0_processed.mat'),...
%                                 'Select Behavior or Recording Session');
% WhereSession = fullfile(SessionPath,WhichSession);

WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q4/Q4_20221109_r0_processed.mat';
WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q9/Q9_20221116_r0_processed.mat';
WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/O3/O3_20210929_r0_processed.mat';
WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q3/Q3_20221019_r0_processed.mat';


%% Load the relevant variables
load(WhereSession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits', ...
               'TimestampAdjust');

%% flag perturbation trials
if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    perturbedTrials = x(~strcmp(TrialInfo.Perturbation(x,1),'OL-Template'));
else
    perturbedTrials = [];
end

%% Change the trial vector to include the OdorON duration before trial start     
for i = 1:size(TrialInfo.TrialID,2)
    OdorTrace = Traces.Trial{i};
    idx2 = find(OdorTrace((SampleRate*startoffset):end)>0,1,'first') + (SampleRate*startoffset) - 1;
    idx1 = idx2 + round(TrialInfo.OdorStart(i,1)*SampleRate);
    
    if idx1<1
        %keyboard;
        
        % for the current trial
        OdorTrace(OdorTrace>0) = 1;
        OdorTrace(1:idx2,1) = 1;
        Traces.Odor{i} = OdorTrace*TrialInfo.Odor(i,1);
        
        % previous trial
        prevOdorTrace = Traces.Odor{i-1};
        prevOdorTrace(end+idx1:end) = TrialInfo.Odor(i,1);
        Traces.Odor{i-1} = prevOdorTrace;
    else
        OdorTrace(OdorTrace>0) = 1;
        OdorTrace(idx1:idx2,1) = 1;
        Traces.Odor{i} = OdorTrace*TrialInfo.Odor(i,1);
    end
    
    % Binarize TrialTrace and Flag Perturbations
    TrialTrace = Traces.Trial{i};
    TrialTrace(TrialTrace>0) = 1;
    
    if ~isempty(find(perturbedTrials==i))
        TrialTrace = TrialTrace*-1;
    end
    Traces.Trial{i} = TrialTrace;
    
end
           
%% Concatenate traces
whichTrials = 1:length(TrialInfo.TrialID);
traceOverlap = SampleRate*startoffset;
whichTraces{1} = 'Lever'; whichTraces{2} = 'Motor'; whichTraces{3} = 'Sniffs';
whichTraces{4} = 'Trial'; whichTraces{5} = 'Odor';
whichTraces{6} = 'Rewards'; whichTraces{7} = 'Timestamps';

for j = 1:size(whichTraces,2)
    temp = cellfun(@(x) ...
        x(1:end-traceOverlap), Traces.(whichTraces{j})(whichTrials), ...
        'UniformOutput', false);
    TracesOut.(whichTraces{j}) = {[cell2mat(temp(:)); ...
        Traces.(whichTraces{j}){whichTrials(end)}(end-traceOverlap+1:end,1)]};
end

% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% snity check for timestamp drops
if any(round(diff(Timestamps),3,'decimal')~=0.002)
    keyboard;
end

% convert the behavior timestamps to OEPS base
TracesOut.Timestamps{1} = Timestamps + TimestampAdjust.ClosedLoop;

%% Get the timestamp for Closed Loop End
SessionLength = ceil(TrialInfo.SessionTimestamps(end,2) + TimestampAdjust.ClosedLoop); % in seconds

%% Organize data for Wolf separately

path = fileparts(fileparts(fileparts(WhereSession)));
if ~exist(fullfile(path,'WDW'),'dir')
    mkdir(fullfile(path,'WDW'));
end
[~,filename] = fileparts(WhereSession);
filename = [filename,'.mat'];
savepath = '/home/priyanka/Desktop';
save(fullfile(savepath,'forWDW',filename),'TracesOut','SingleUnits');
