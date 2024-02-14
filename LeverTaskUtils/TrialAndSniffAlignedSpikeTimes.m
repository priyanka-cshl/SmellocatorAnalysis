function [TrialAlignedSniffs, SniffAlignedSpikes, TrialAlignedSpikes, TetrodeOrder, Events, Events_phase, TrialInfo] = TrialAndSniffAlignedSpikeTimes(SingleUnits,TTLs,TrialInfo,PerturbationTag)
 
N = size(SingleUnits,2); % total units

if nargin < 4
    PerturbationEvents = [];
else
    % to get comparable alignment time-points as perturbartion start in ctrl trials
    [PerturbationEvents] = FindPutativePerturbationStart(TrialInfo, PerturbationTag);
    % for halts:    [nx3]:[HaltStart HaltEnd HaltLocation] 
    % for offsets:  [nx3]:[[OffsetStart FeedbackResume OffsetLocation]
    % for reverals: [nx1]:[[1 if reversed, 0 otherwise]
    % n = trials
end

%% Sniff alignment
% get trial start and stop times
trialtimes = TrialInfo.SessionTimestamps(:,[1 2 1]);
% recompute Trialstart in behavior timebase - only need this because of
% sample drops at trial start
funkytrials = find(TrialInfo.TimeStampsDropped);
if ~isempty(funkytrials)
    trialtimes(funkytrials,1) = trialtimes(funkytrials,2) - TTLs.Trial(funkytrials,3);
end
trialtimes(:,3) = [];

% the sniff time-points (same time base as the trial times)
load(TrialInfo.SessionPath,'SniffTS'); % sniff timestamps in behavior timebase

TrialAlignedSniffs = GetTrialAlignedSniffs(trialtimes,SniffTS);
% type help GetTrialAlignedSniffs to check notes on the output

%% Align Spikes to Trials and also parse by sniffs
sniffwarpmethod = 3; % fixed latency
nTrials     = size(TrialInfo.TrialID,2); 
for whichunit = 1:N % every unit
    thisUnitspikes = SingleUnits(whichunit).trialalignedspikes;
    trialtags = SingleUnits(whichunit).trialtags;
    TetrodeOrder(whichunit,:) = [SingleUnits(whichunit).tetrode SingleUnits(whichunit).id]; % tetrode and phy cluster ID

    for whichtrial = 1: nTrials % every trial
        thisTrialspikes = thisUnitspikes(trialtags==whichtrial);
        if whichtrial>1
            previousTrialspikes = thisUnitspikes(trialtags==whichtrial-1) - ...
                (TTLs.Trial(whichtrial,1) - TTLs.Trial(whichtrial-1,1));
            
            % ignore spikes that happened before previous trial off
            prevToff = TTLs.Trial(whichtrial-1,2) - TTLs.Trial(whichtrial,1);
            if ~isempty(previousTrialspikes)
                previousTrialspikes(previousTrialspikes<prevToff) = [];
            end
        else
            previousTrialspikes = thisUnitspikes(trialtags==whichtrial-1) - ...
                (TTLs.Trial(whichtrial,1) - 0);
        end
                
        % convert every spike to sniff phase
        myspikes = horzcat(previousTrialspikes',thisTrialspikes');
        TrialAlignedSpikes{whichtrial,whichunit} = { myspikes };
        SniffAlignedSpikes{whichtrial,whichunit} = { WhichSniffPhase(myspikes,TrialAlignedSniffs{whichtrial},'warpmethod',sniffwarpmethod) };
        
        if whichunit == 1
            x1 = TTLs.Trial(whichtrial,1);
            x2 = TTLs.Trial(whichtrial,2) + 0.2; % next trial cannot start sooner than that
            reward = TTLs.Reward(find((TTLs.Reward(:,1)>x1)&(TTLs.Reward(:,1)<x2)),1);
            if isempty(reward)
                reward = NaN;
            end
            % OdorStart, Reward, TrialOFF, First TZ entry
            Events(whichtrial,1:4) = [TTLs.Trial(whichtrial,4), ...
                (reward(1) - x1), ...
                TTLs.Trial(whichtrial,3), ...
                TrialInfo.TargetEntry(whichtrial,1)];
            
            if ~isempty(PerturbationEvents)
                Events(whichtrial,5) = PerturbationEvents(whichtrial,1)/SampleRate;
            end
            
            % alignEventTimes to Sniffs
            Events_phase(whichtrial,:) = WhichSniffPhase(Events(whichtrial,:),TrialAlignedSniffs{whichtrial},'warpmethod',sniffwarpmethod);
            
            % also recalculate InZone periods in sniff base
            if ~isempty(TrialInfo.InZone{whichtrial})
                TrialInfo.InZonePhase{whichtrial}(:,1) = WhichSniffPhase(TrialInfo.InZone{whichtrial}(:,1),TrialAlignedSniffs{whichtrial});
                TrialInfo.InZonePhase{whichtrial}(:,2) = WhichSniffPhase(TrialInfo.InZone{whichtrial}(:,2),TrialAlignedSniffs{whichtrial});
            end

        end
        
    end
end

end