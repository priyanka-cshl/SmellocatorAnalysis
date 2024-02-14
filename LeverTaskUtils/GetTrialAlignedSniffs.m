function [TrialAlignedSniffs] = GetTrialAlignedSniffs(trialtimes, SniffTS)
% helper function to express Sniff timestamps relative to trial start
% Input: trialtimes [n x 4(2)]:
%           [trialstart trialend prevtrialend nexttrialstart]
%      : SniffTS [n x 3]
%           [InhalationStart InhalationEnd NextInhalationStart]
% Output: {t x 1} where each {} = [nx14]; t = trials, 
%   n = sniffs within each trial (from previous trial off to this trial end)
%   Cols [1 2 3 6], [7 8 9 10], [11 12 13 14] are 
%   InhStart InhEnd InhStart OdorLocation for current, prev and next sniff
%   Col 4 = sniff index w.r.t. trial start, where 1 = first sniff
%   Col 5 = sniff type (from -1 to 2) for current sniff
%       -1 (ITI)
%        0 (in trial, before first significant TZ entry)
%        1 (in trial, after first significant TZ entry)
%        2 (in target zone)

nTrials = size(trialtimes,1); % total no. of trials
if size(trialtimes,2) == 2
    % this will happen when passing trialtimes from the closed loop session
    % compute the previous trial off and next trial ON times
    trialtimes(1,3)         = 0; % session start
    trialtimes(2:end,3)     = trialtimes(1:end-1,2);
    trialtimes(1:end-1,4)   = trialtimes(2:end,1);
    trialtimes(end,4)       = trialtimes(end,2) + 1; % 1 sec after the last trial
end

if ~isempty(SniffTS)
    % get a cell array trial aligned sniffs
    for whichtrial = 1: nTrials % every trial
%         if whichtrial==1
%             first_inhalation = 2;
%         else
%             first_inhalation = find(SniffTS(:,1)>=trialtimes(whichtrial-1,2),1,'first'); % after prev trial OFF
%         end
%         last_inhalation  = find(SniffTS(:,2)<=trialtimes(whichtrial,2),1,'last'); % before trial end
        
        first_inhalation = find(SniffTS(:,1)>=trialtimes(whichtrial,3),1,'first'); % after prev trial OFF
        first_inhalation = max(first_inhalation,2); % need the minimum to be 2, such that we can compute the previous sniff TS
        last_inhalation  = find(SniffTS(:,2)<=trialtimes(whichtrial,2),1,'last'); % before trial end
        
        mysniffs = SniffTS(first_inhalation:last_inhalation,1:3) - trialtimes(whichtrial,1); % -ve timestamps are in ITI
        mysniffs(:,4) = (1:size(mysniffs,1))' - numel(find(mysniffs(:,1)<=0));
        mysniffs(:,6) = SniffTS(first_inhalation:last_inhalation,4); % odor location 
        
        % flag sniffs where entire inhalation period was in the target
        % zone?
        
        % old way - use TargetZone entry times
%         if ~isempty(TrialInfo.InZone{whichtrial})
%             whichsniffs = find(mysniffs(:,1)>=TrialInfo.TargetEntry(whichtrial,1));
%             mysniffs(whichsniffs,5) = 1;
%             for entries = 1:size(TrialInfo.InZone{whichtrial},1)
%                 stay = TrialInfo.InZone{whichtrial}(entries,:);
%                 whichsniffs = intersect( (find(mysniffs(:,1)>=stay(1))) , ...
%                     (find(mysniffs(:,2)<=stay(2))) );
%                 mysniffs(whichsniffs,5) = 2;
%             end
%         end
        
        % better way - use the location
        whichsniffs = intersect(find(mysniffs(:,3)>0), find(abs(mysniffs(:,6))<=9));
        
        mysniffs(mysniffs(:,3)<=0,5) = -1;
        
        % also note the sniff before and sniff after
        mysniffs = horzcat(mysniffs, ...
            SniffTS(first_inhalation-1:last_inhalation-1,1:3) - trialtimes(whichtrial,1), ...
            SniffTS(first_inhalation-1:last_inhalation-1,4), ...
            SniffTS(first_inhalation+1:last_inhalation+1,1:3) - trialtimes(whichtrial,1), ...
            SniffTS(first_inhalation+1:last_inhalation+1,4) );
        
        TrialAlignedSniffs{whichtrial} = mysniffs;
        
    end
end

end