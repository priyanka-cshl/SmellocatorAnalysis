function [ReplayTraces] = FixFunkyReplays(ReplayTraces,ReplayLengths,templateindex)
%% function to fix the issue with some AON open loop sessions     
% where some replays have corrupt and extra samples
%%

if nargin<3
    templateindex = 1;
end

global SampleRate; % = 500; % samples/second
global startoffset; % = 1; % seconds
traceOverlap = SampleRate*startoffset + 1; 

% step 1: find the culprit replays
buggy_replays   = find(abs(ReplayLengths' - median(ReplayLengths))>5);
[~,OK_replay]   = min(ReplayLengths' - median(ReplayLengths));
OK_replay = ReplayTraces.Motor{templateindex}(:,OK_replay); % backwards from trial end
for b = 1:numel(buggy_replays)
    extra_samples = ReplayLengths(buggy_replays(b)) - median(ReplayLengths);
    if extra_samples>0
        % buggy_replay have extra samples and possibly some
        % corrupt ones surrounding it
        buggy_backwards = ReplayTraces.Motor{templateindex}(:,buggy_replays(b));
        % traces will initially align (from trial end) and then
        % error will increase because of corrupt samples in the
        % buggy replay - find the inflexion point
        smoothing_window = 25; % 50 ms window for running average
        mismatch_point = find(abs(smooth(buggy_backwards(traceOverlap:end) - OK_replay(traceOverlap:end),25))>2,1,'first');
        
        
        
        % try to align the remaining trace by ignore some samples  
        temp1  = OK_replay(mismatch_point:(end-traceOverlap));
        temp2  = buggy_backwards((mismatch_point+extra_samples):(end-traceOverlap));
        [r,lags] = xcorr(temp1,temp2,4*extra_samples);
        if numel(find(r==max(r))) == 1 % unique maximum correlation
            samps_to_discard = traceOverlap - 1 + mismatch_point + (1:(extra_samples - lags(find(r==max(r)))));
            adjusted_vec = buggy_backwards;
            adjusted_vec(samps_to_discard,:) = [];
            % sanity check
            min_length = min(length(adjusted_vec),length(OK_replay));
            match_quality = corrcoef([adjusted_vec(1:min_length,1) OK_replay(1:min_length,1)]);
            if match_quality(1,2)>=0.97
                % do the corrupt samples fall in the ITI period?
                ReplayTraces.Timestamps{templateindex}(samps_to_discard,buggy_replays(b)) = ...
                    -ReplayTraces.Timestamps{templateindex}(samps_to_discard,buggy_replays(b));
                disp(['fixed buggy replay #', num2str(b)]);
                
            else
                disp('having trouble fixing buggy replays');
                keyboard;
            end
        end
        
    else
        disp('having trouble fixing buggy replays');
        keyboard;
    end
end

end