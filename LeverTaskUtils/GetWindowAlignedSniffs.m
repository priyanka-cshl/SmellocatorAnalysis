function [TrialAlignedSniffs, RefTimeStamps] = GetWindowAlignedSniffs(trialtimes, SniffTS, trialinfo)
% helper function to express Sniff timestamps relative to trial/event start
% Input: trialtimes [n x 4(2)]:
%           [trialstart trialend prevtrialend nexttrialstart]
%      : SniffTS [n x 3]
%           [InhalationStart InhalationEnd NextInhalationStart]
% Output: TrialAlignedSniffs: {t x 1} where, each {} = [nx16]; 
%           t = trials, n = sniffindex within a trial
%           Cols 1-4 are for details of the current sniff
%           Col 1: TrialID
%           Col 2: Odor/Air State: -1 (Air OFF), 0-3 (air, odor 1,2,3)
%           Col 3: sniff index w.r.t. event start, where 1 = first sniff
%           Col 4: sniff type (from -1 to 2) for current sniff
%                   -1 (ITI)
%                    0 (in trial, before first significant TZ entry)
%                    1 (in trial, after first significant TZ entry)
%                    2 (in target zone)
%           Cols 5-11: are timestamps for prev, current and next sniffs
%               PrevInhStart PrevInhEnd 
%               CurrentInhStart CurrentInhEnd
%               NextInhStart NextInhEnd NextExhEnd
%           Cols 12-14: MotorLocation for prev, current, and next sniff
%       : RefTimeStamp: [t x 1], where each entry is the reference window
%           start time that was subtracted from the raw sniff timestamps


nTrials = size(trialtimes,1); % total no. of trials
nSniffs = size(SniffTS,1); % total no. of sniffs

if ~isempty(SniffTS)
    % get a cell array trial aligned sniffs
    for whichtrial = 1: nTrials % every trial
        
        first_inhalation = find(SniffTS(:,1)>=trialtimes(whichtrial,3),1,'first'); % after prev trial OFF
        first_inhalation = max(first_inhalation,2); % need the minimum to be 2, such that we can compute the previous sniff TS
        
        last_inhalation  = find(SniffTS(:,2)<=trialtimes(whichtrial,2),1,'last'); % before trial end
        
        % recmpoute snif timestamps w.r.t. Trial/Window Start 
        % -ve timestamps are in ITI (Pre-Event)
        currentsniff        = SniffTS(first_inhalation:last_inhalation,1:3)     - trialtimes(whichtrial,1);
        previoussniff       = SniffTS(first_inhalation-1:last_inhalation-1,1:3) - trialtimes(whichtrial,1); 
        nextsniff           = SniffTS(first_inhalation+1:last_inhalation+1,1:3) - trialtimes(whichtrial,1);
        
        snifftimestamps     = [previoussniff(:,1:2) currentsniff nextsniff(:,2:3)];
        
        motorlocations      = [SniffTS(first_inhalation-1:last_inhalation-1,4), ...
                               SniffTS(first_inhalation:last_inhalation,4),     ...
                               SniffTS(first_inhalation+1:last_inhalation+1,4)  ];
                           
        % details about the current sniff
        sniffdetails        = zeros(size(motorlocations,1),4);
        sniffdetails(:,1)   = whichtrial; % trial ID
        sniffdetails(:,2)   = 0;  % Odor State
        if isfield(trialinfo,'TuningTrialID') % sniff from passive tuning trials
            sniffdetails(:,2) = trialinfo.Odor(whichtrial,1); % odor state has already been computed in PreprocessSpikesAndSniffs.m
        elseif isfield(trialinfo,'OdorStart') % close loop trials or replays
            odorON = trialinfo.OdorStart(whichtrial,1); % seconds w.r.t. trial start
            sniffdetails((currentsniff(:,1)< odorON),2) = -1; % ITI (air OFF)
            sniffdetails((currentsniff(:,1)>=odorON),2) = trialinfo.Odor(whichtrial);
        end
        % sniff index (w.r.t. trial start)
        sniffdetails(:,3) = (1:size(sniffdetails,1))' - numel(find(currentsniff(:,1)<=0));
         
        % sniff type
        sniffdetails(:,4)   = 0; 
        % ITI sniffs
        if isfield(trialinfo,'OdorStart')
            odorON = trialinfo.OdorStart(whichtrial,1); % seconds w.r.t. trial start
            sniffdetails((currentsniff(:,1)< odorON),4) = -1; % ITI (air OFF)
        else
            sniffdetails(currentsniff(:,1)<0,4) = -1;
        end
        if isfield(trialinfo,'TargetEntry')
            sniffdetails(currentsniff(:,1)>trialinfo.TargetEntry(whichtrial),4) = 1; % after Target Entry
        end
        % flag sniffs where entire inhalation period was in the target zone
        % use the motor location
        whichsniffs = intersect(find(currentsniff(:,1)>=0), find(abs(motorlocations(:,2))<=9));
        sniffdetails(whichsniffs,4) = 2;       
        
        TrialAlignedSniffs{whichtrial}  = [sniffdetails snifftimestamps motorlocations];
        RefTimeStamps(whichtrial,1)     = trialtimes(whichtrial,1); 
    end
end

end