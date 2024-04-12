function [TuningInfo] = ParseTuningSubtrials(TuningTTLs, TuningSequence, TuningSettings)

% Extract timestamps for periods with odor On, Air On and Air Off
% unlike for behavior trials - no bufferWindows for before odor etc
% since the pre/post periods differ a lot in tuning settings

% 1: determine passive tuning type
if size(TuningSequence,1)>2 && size(TuningSequence,2) == 2
    tuningtype = 0; % old style passive tuning
elseif any(TuningSequence(:,1)==800)
    tuningtype = 1; % pseudorandom tuning
end

% 2: determine if Air was ON or OFF during ITI
if length(TuningSettings) == 12
    ITIAirState = TuningSettings(12);
else
    ITIAirState = 0;
end

% TuningSettings = Timestamps w.r.t. trial start ...
% [~ ~  motor-settle pre-odor odor purge post-odor ITI]
TuningSettings(3:8) = TuningSettings(3:8)/1000; % convert to seconds
if tuningtype % pseudorandom tuning
    % get general parameters for each pesudorandom trial
    % ie. TS w.r.t trial start at which location updates happen
    
    % due to a bug, for the first location, odor turns on actually after
    % the location update, so need to keep track of that extra transition
    % also, if air is ON during ITI, pre-trial period can be included in
    % the baseline for location #1, before odor actually turns ON
    
    % total locations covered
    nLocations     = size(TuningSequence,2) - 2;
    StateShifts = zeros(nLocations+1,4); % [tstart tend locationID odorstate]
    
    % two transitions to account for before the update to location #2
    
    % first transition is until odor actually turns on
    if ITIAirState
        StateShifts(1,1) = -TuningSettings(1,3); % w.r.t. trial start (settle)
    else
        StateShifts(1,1) = 0; % trial start (before that Air was OFF)
    end
    StateShifts(1,2) = TuningSettings(1,4); % w.r.t. trial start (pre-odor)
    StateShifts(1,3) = 1; % first location from the location vector for that trial
    StateShifts(1,4) = 0; % Air is ON, but not odor
    
    % second transition (same location as above, but now odor is on)
    StateShifts(2,1) = TuningSettings(1,4); % w.r.t. trial start (pre-odor)
    StateShifts(2,2) = StateShifts(2,1) + TuningSettings(1,5); % add odor
    StateShifts(2,3) = 1; % first location from the location vector that trial
    StateShifts(2,4) = 1; % Odor is ON
    
    % now update for all remaining locations
    
    Loc_duration    = sum(TuningSettings(1,[3,5])); % settle + odor
    StateShifts(3:end,2) = StateShifts(2,2) + cumsum(repmat(Loc_duration,nLocations-1,1));
    StateShifts(3:end,1) = StateShifts(3:end,2) - Loc_duration;
    StateShifts(3:end,3) = 2:nLocations;
    StateShifts(3:end,4) = 1; % odor is ON
    
    % add one or two more transitions for the post-odor period
    nTransitions = size(StateShifts,1);
    StateShifts(nTransitions+1,3) = nLocations; % last location from the location vector for that trial
    StateShifts(nTransitions+1,4) = 0; % Air is ON, but not odor
    StateShifts(nTransitions+1,1) = StateShifts(nTransitions,2) ; % start of post-odor period
    if ITIAirState
        StateShifts(nTransitions+1,2) = StateShifts(nTransitions+1,1) + ...
            sum(TuningSettings(1,7:8)); % post-odor + ITI (both have air On, odor OFF and same location
    else
        StateShifts(nTransitions+1,2) = StateShifts(nTransitions+1,1) + ...
            TuningSettings(1,7); % just the post-odor, after that air turns OFF
        
        % add another transition for the ITI period
        StateShifts(nTransitions+2,1) = StateShifts(nTransitions+1,2) ; % start of ITI
        StateShifts(nTransitions+2,2) = StateShifts(nTransitions+2,1) + ...
            TuningSettings(1,8); % ITI only
        StateShifts(nTransitions+2,3) = nLocations; % last location from the location vector for that trial
        StateShifts(nTransitions+2,4) = -1; % Air is OFF
        
        % Also add another transition in the very beginning
        % before trial start - when air was OFF
        % [-settle trialstart location#1 AirOFF]
        StateShifts = vertcat([-TuningSettings(1,3) 0 1 -1], StateShifts);
        
    end
else
    % regular tuning
    % get general parameters for each tuning trial
    % ie. TS w.r.t trial start at which odor/air transitions happen
    
    nStates = 3 + 2*~ITIAirState;
    % [settle PreOdor Odor PostOdor ITI]
    % if air is ON during ITI, settle = PreOdor = Postodor = ITI in terms of location and valve states
    
    StateShifts = zeros(3,4); % [tstart tend locationID odorstate]
    
    % first transition is until odor actually turns on
    if ITIAirState
        StateShifts(1,1) = -TuningSettings(1,3); % w.r.t. trial start (settle)
    else
        StateShifts(1,1) = 0; % trial start (before that Air was OFF)
    end
    StateShifts(1,2) = TuningSettings(1,4); % w.r.t. trial start (pre-odor)
    StateShifts(1,3) = 1; % first location from the location vector for that trial
    StateShifts(1,4) = 0; % Air is ON, but not odor
    
    % second transition is to post-odor
    StateShifts(2,1:2) = StateShifts(1,2) + [0 TuningSettings(1,5)]; % odor
    StateShifts(2,3:4) = [1 1]; % same location, odor ON
    
    % third transition is to ITI
    if ITIAirState
        StateShifts(3,1:2) = StateShifts(2,2) + [0 sum(TuningSettings(1,6:8))]; % purge + postodor + ITI
    else
        StateShifts(3,1:2) = StateShifts(2,2) + [0 sum(TuningSettings(1,6:7))]; % purge + postodor
    end
    StateShifts(3,3:4) = [1 0]; % same location, odor OFF, air ON
    
    if ~ITIAirState
        StateShifts(4,1:2) = StateShifts(3,2) + [0 sum(TuningSettings(1,8))]; % ITI
        StateShifts(4,3:4) = [1 -1]; % same location, Air OFF
        
        StateShifts(5,1)   = -TuningSettings(1,3); % w.r.t. trial start (settle)
        StateShifts(5,2)   = 0; % trialstart
        StateShifts(5,3:4) = [1 -1]; % same location, odor OFF, air ON
        
        StateShifts = circshift(StateShifts,1);
    end
end

%% parse to subtrials
% ignore any non-tuning trials
if tuningtype
    IgnoreTrials = find(TuningSequence(:,1) ~= 800);
else
    IgnoreTrials = find(TuningSequence(:,1) > 800);
end
% append an extra column for keeping track of
% trial ID (w.r.t. to Tuning TTLs)
TuningTTLs(:,end+1) = 1:size(TuningTTLs,1);
TuningTTLs(IgnoreTrials,:) = [];
TuningSequence(IgnoreTrials,:) = [];
nTrials = size(TuningTTLs,1); % no. of valid pseudorandom tuning trials

AllSubtrials = [];
for N = 1:nTrials % every tuning trial
    thisTrialStart = TuningTTLs(N,1); % trial start in OEPS base
    thisTrialTimes = StateShifts(:,1:2) + thisTrialStart; % tstart and tstop
    thisTrialTimes(:,3:4) = thisTrialTimes(:,1:2); % no buffer window for prev trial off and next trial start
    if ~tuningtype
        thisTrialTimes(:,5) = TuningTTLs(N,7)*StateShifts(:,3); % Motor Locations
        thisTrialTimes(:,6) = TuningTTLs(N,5)*StateShifts(:,4); % odor IDs
    else
        thisTrialTimes(:,5) = TuningSequence(N, 2 + (StateShifts(:,3))); % Motor Locations
        thisTrialTimes(:,6) = TuningSequence(N,2)*StateShifts(:,4); % odor IDs
    end
    thisTrialTimes(:,6) = max(thisTrialTimes(:,6),-1);
    thisTrialTimes(thisTrialTimes(:,6)>0,6) = thisTrialTimes(thisTrialTimes(:,6)>0,6) - 1; % change odor indices from 1-4 to 0-3
    thisTrialTimes(:,7) = TuningTTLs(N,8); % trialID w.r.t. TTLs.Trial
    thisTrialTimes(:,8) = TuningTTLs(N,12) + (1:size(StateShifts,1))/100; % trialID w.r.t. TuningTTLs
    AllSubtrials = vertcat(AllSubtrials, thisTrialTimes);
end

% Useful stuff for parsing subtrials
nSubTrials = size(AllSubtrials,1);
TuningInfo.trialtimes       = AllSubtrials(:,1:4);
TuningInfo.Location         = AllSubtrials(:,5);
TuningInfo.Odor             = AllSubtrials(:,6);
TuningInfo.Perturbed        = zeros(nSubTrials,1);
TuningInfo.Perturbation     = NaN;
TuningInfo.TrialID          = (1:nSubTrials)';
TuningInfo.OriginalTrialID  = AllSubtrials(:,7);
TuningInfo.TuningTrialID    = AllSubtrials(:,8);

end


