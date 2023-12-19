function [] = FixFunkyReplays(MotorTrace,RewardTrace,ReplayLengths)
%% function to fix the issue with some AON open loop sessions     
% where some replays have corrupt and extra samples
%%

global SampleRate; % = 500; % samples/second
global startoffset; % = 1; % seconds
traceOverlap = SampleRate*startoffset + 1; 

% step 1: find the culprit replays
buggy_replays   = find(abs(ReplayLengths' - median(ReplayLengths))>5);
[~,OK_replay]   = min(ReplayLengths' - median(ReplayLengths));
OK_replay = MotorTrace(:,OK_replay); % backwards from trial end
for b = 1:numel(buggy_replays)
    extra_samples = ReplayLengths(buggy_replays(b)) - median(ReplayLengths);
    if extra_samples>0
        % buggy_replay have extra samples and possibly some
        % corrupt ones surrounding it
        buggy_backwards = MotorTrace(:,buggy_replays(b));
        % traces will initially align (from trial end) and then
        % error will increase because of corrupt samples in the
        % buggy replay - find the inflexion point
        smoothing_window = 25; % 50 ms window for running average
        mismatch_point = find(abs(smooth(buggy_backwards(traceOverlap:end) - OK_replay(traceOverlap:end),25))>2,1,'first');
        
        
        
        % try to align the remaining trace by ignore some samples  
        temp1  = OK_replay((mismatch_point+extra_samples):(end-traceOverlap));
        temp2  = buggy_backwards(mismatch_point:end);
        [r,lags] = xcorr(temp1,temp2,2*extra_samples);
        
        
    else
        keyboard;
    end
end

end