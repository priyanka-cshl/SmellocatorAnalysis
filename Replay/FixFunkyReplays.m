function [ReplayTraces] = FixFunkyReplays(ReplayTraces,ReplayLengths,ReplayTTLs,templateindex)
%% function to fix the issue with some AON open loop sessions     
% where some replays have corrupt and extra samples
%%

if nargin<4
    templateindex = 1;
end

global SampleRate; % = 500; % samples/second
global startoffset; % = 1; % seconds
traceOverlap = SampleRate*startoffset + 1; 
whichTraces = fieldnames(ReplayTraces);

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
        mismatch_point = find(abs(smooth(buggy_backwards(traceOverlap:end) - OK_replay(traceOverlap:end),smoothing_window))>2,1,'first');
        
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
            if match_quality(1,2)>=0.95
                disp(['fixed buggy replay #', num2str(buggy_replays(b))]);
                % do the corrupt samples fall in the ITI period?
                % remap timestamps relative to replay end
                thisReplayTTLs = ReplayTTLs.OdorValve{buggy_replays(b)};
                thisReplayEnd = ReplayTraces.Timestamps{templateindex}(traceOverlap,buggy_replays(b));
                tempTS = ReplayTraces.Timestamps{templateindex}(:,buggy_replays(b)) - thisReplayEnd + thisReplayTTLs(end,2);
                bug_ended = tempTS(samps_to_discard(1),1);
                bug_began = tempTS(samps_to_discard(end),1);
                afterTrial = find(thisReplayTTLs(:,2)>=bug_began, 1, 'first');
                if find(thisReplayTTLs(:,1)>=bug_ended, 1, 'first') == afterTrial
                    disp(['buggy samples were in the ITI after subtrial# ',num2str(afterTrial-1)]);
                    ReplayTraces.Corrupt{templateindex}(buggy_replays(b),:) = [afterTrial-1 0 bug_began bug_ended numel(samps_to_discard)];
                else
                    disp(['buggy samples were within subtrial# ',num2str(afterTrial)]);
                    ReplayTraces.Corrupt{templateindex}(buggy_replays(b),:) = [afterTrial 1 bug_began bug_ended numel(samps_to_discard)];
                end
                
                % edit the Replay traces to chop out the corrupt samples 
                % and just put them in the end 
                % traces will get flipped back so these will actually go in
                % the very beginning of the corrected traces
                for j = 2:size(whichTraces,1) % first is trial IDs, last is corruption
                    temptrace = ReplayTraces.(whichTraces{j}){templateindex}(:,buggy_replays(b));
                    ReplayTraces.(whichTraces{j}){templateindex}(:,buggy_replays(b)) = ...
                            vertcat(temptrace(setdiff(1:length(temptrace),samps_to_discard),:), ...
                                temptrace(samps_to_discard,:));
                            
                    if strcmp(whichTraces{j},'Timestamps')
                        temptrace(samps_to_discard,:) = -temptrace(samps_to_discard,:);
                    end
                    
                    ReplayTraces.Uncorrected.(whichTraces{j}){templateindex}(:,buggy_replays(b)) = temptrace;
                end
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

% extra step to remove trailing samples at the end of each trace (this will
% go to the beginning of the traces when flipped back and cause
% misalignment)
trailing_samps = max(ReplayTraces.Corrupt{templateindex}(:,5)) - 1;
for j = 2:size(whichTraces,1)
    ReplayTraces.(whichTraces{j}){templateindex}(end-trailing_samps:end,:) = [];
end

end