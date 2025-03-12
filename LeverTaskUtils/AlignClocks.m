function [TimestampAdjust] = AlignClocks(WhereSession)
% Calculate the exact time adjustment between behavior and ephys timestamps

% Load relevant stuff
load(WhereSession, 'Traces', 'TrialInfo', ...
    'startoffset', 'SampleRate', ...
    'TTLs', ...
    'TimestampAdjust');

if exist('TimestampAdjust','var')
    if length(TimestampAdjust.ClosedLoop) == 2
        return;
    end
end

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

% calculate time adjustment
ephys = TTLs.Trial(1:numel(TrialInfo.Odor),2);
behavior = Timestamps(find((diff(abs(TracesOut.Trial{1})))==-1))';
myfit = fit(behavior,ephys-behavior,'poly1');

TimestampAdjust.ClosedLoop(2) = myfit.p2; 
TimestampAdjust.ClosedLoop(1) = myfit.p1;
