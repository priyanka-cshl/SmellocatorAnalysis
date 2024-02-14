function [PerturbationEvents] = FindPutativePerturbationStart(TrialInfo, PerturbationTag)
%   helper function to find the time-point (w.r.t. trial start)
%   that is analogous to perturbation start in CONTROL trials
%   for perturbation types - halts, offsets and Rule-reversals

load(TrialInfo.SessionPath,'Traces','SampleRate'); % Traces is used only for halts, SampleRate is needed for all
PerturbationEvents = [];
Threshold   = 4; % Lever threshold at which halt perturbation starts
nTrials     = size(TrialInfo.TrialID,2); % total no. of trials

if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1))) % are there any perturbation trials in this session?
    
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1))); % find trial # for perturbation trials
    
    if regexpi(PerturbationTag, 'Reversal')
        PerturbationEvents(1:nTrials,1) = 0;
        % mark all flipped trials as one
        PerturbationEvents(find(strcmp(TrialInfo.Perturbation(:,1),'RuleReversal')),1) = 1;
        return;
    end
    
    PerturbationEvents = repmat([NaN NaN NaN],nTrials,1);
    
    PerturbedTrials = [];
    
    if regexpi(PerturbationTag, 'halt')
        %% Halt-flip trials
        PerturbedTrials = x(find(cellfun(@(x) ~isempty(regexpi(x,'halt')), TrialInfo.Perturbation(x,1))));
        
        PerturbationParams = TrialInfo.Perturbation{PerturbedTrials(1),2}; % start and stop Idx w.r.t. TrialStart, and halted odor location
        Duration = PerturbationParams(2) - PerturbationParams(1);
    end
    
    if regexpi(PerturbationTag, 'offset')
        %% Offset trials
        PerturbedTrials = x(find(cellfun(@(x) ~isempty(regexpi(x,'Offset-II')), TrialInfo.Perturbation(x,1))));
        
        PerturbationParams = TrialInfo.Perturbation{PerturbedTrials(1),2}; % offset and perturbation start w.r.t. TrialStart, and offset odor location
        Duration = PerturbationParams(2) - PerturbationParams(1);
    end
    
    % for all control trials, find perturbation start point
    for i = 1:nTrials
        
        if ~ismember(i,PerturbedTrials) % control trials
            
            if regexpi(PerturbationTag, 'halt')
                % find threshold crossing - when would a halt get triggered
                f = find(Traces.Lever{i}(TrialInfo.TimeIndices(i,1):end)<Threshold,1,'first');
                if ~isempty(f)
                    f = f + TrialInfo.TimeIndices(i,1) - 1;
                    if abs(Traces.Lever{i}(f-1)-4)<abs(Traces.Lever{i}(f)-4)
                        f = f-1;
                    end
                    f = f - TrialInfo.TimeIndices(i,1);
                    PerturbationEvents(i,:) = [f f+Duration NaN];
                end
            end
            
            if regexpi(PerturbationTag, 'offset')
                % get all target zone stays
                thisTrialHolds = TrialInfo.HoldSettings(i,2:3); % contiguous, aggegate in seconds
                AllHolds = diff(TrialInfo.InZone{i}',1);
                                
                if ~isempty(AllHolds)
                    whichSegment = find(AllHolds>=0.9*thisTrialHolds(1),1,'first'); % any contiguous hold that would trigger offset perturbation
                    if isempty(whichSegment)
                        whichSegment = ...
                            find(cumsum(AllHolds)>=0.9*thisTrialHolds(2),1,'first'); % any aggregate hold that would trigger offset perturbation
                        if ~isempty(whichSegment)
                            deltaT = sum(AllHolds)-0.9*thisTrialHolds(2);
                            PerturbationEvents(i,1) = SampleRate* ...
                                (TrialInfo.InZone{i}(whichSegment,2) - deltaT); % in indices w.r.t. TrialStart
                        end
                    else
                        PerturbationEvents(i,1) = SampleRate* ...
                            (TrialInfo.InZone{i}(whichSegment,1) + thisTrialHolds(1)*0.9); % in indices w.r.t. TrialStart
                    end
                end
            end
            
        else % perturbed trial
            if ~isempty(TrialInfo.Perturbation{i,2})
                PerturbationEvents(i,:) = TrialInfo.Perturbation{i,2};
            end
        end
        
    end
end

end