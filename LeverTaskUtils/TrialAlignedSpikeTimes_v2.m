function [AlignedSpikes, Events, TetrodeOrder, SniffAlignedSpikes, Events_phase, TrialInfo] = ...
    TrialAlignedSpikeTimes_v2(SingleUnits,TTLs,nTrials,TrialInfo,MySession,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('sniffwarpmethod', 0, @(x) isnumeric(x)); % 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency

% extract values from the inputParser
params.parse(varargin{:});
sniffwarpmethod = params.Results.sniffwarpmethod;

N = size(SingleUnits,2); % total units

PerturbationEvents = [];
if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1))); % find perturbation trials
    
    if any(strcmp(TrialInfo.Perturbation(x,1),'Halt-Flip')) || ...
            any(strcmp(TrialInfo.Perturbation(x,1),'Halt-Flip-Template')) 
        
        HaltTrials = [find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')); ...
            find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip-Template')) ]; % halt flip trials
        
        load(MySession,'Traces','SampleRate');
        PerturbationParams = TrialInfo.Perturbation{HaltTrials(1),2}; % start and stop Idx w.r.t. TrialStart, and halted odor location
        Duration = PerturbationParams(2) - PerturbationParams(1);
        Threshold = 4; % Lever threshold at which perturbation starts
        
        % for all control trials, find perturbation start point
        for i = 1:nTrials
            if ~ismember(i,HaltTrials)
                f = find(Traces.Lever{i}(TrialInfo.TimeIndices(i,1):end)<Threshold,1,'first');
                if ~isempty(f)
                    f = f + TrialInfo.TimeIndices(i,1) - 1;
                    if abs(Traces.Lever{i}(f-1)-4)<abs(Traces.Lever{i}(f)-4)
                        f = f-1;
                    end
                    f = f - TrialInfo.TimeIndices(i,1);
                    PerturbationEvents(i,:) = [f f+Duration NaN];
                else
                    PerturbationEvents(i,:) = [NaN NaN NaN];
                end
            else
                if ~isempty(TrialInfo.Perturbation{i,2})
                    PerturbationEvents(i,:) = TrialInfo.Perturbation{i,2};
                else
                    PerturbationEvents(i,:) = [NaN NaN NaN];
                end
            end
        end
    end
    
    if any(strcmp(TrialInfo.Perturbation(x,1),'Offset-II')) || ...
            any(strcmp(TrialInfo.Perturbation(x,1),'Offset-II-Template'))
        
        OffsetTrials = [find(strcmp(TrialInfo.Perturbation(:,1),'Offset-II')); ...
            find(strcmp(TrialInfo.Perturbation(:,1),'Offset-II-Template')) ];
        
        load(MySession,'SampleRate');
        PerturbationParams = TrialInfo.Perturbation{OffsetTrials(1),2}; % offset and perturbation start w.r.t. TrialStart, and offset odor location
        Duration = PerturbationParams(2) - PerturbationParams(1);
        
        % for all control trials, find perturbation start point
        % time-point when 0.9x of the hold times are reached
        for i = 1:nTrials
            if ~ismember(i,OffsetTrials)
                
                % get all target zone stays
                thisTrialHolds = TrialInfo.HoldSettings(i,2:3); % contiguous, aggegate in seconds
                AllHolds = diff(TrialInfo.InZone{i}',1);
                
                PertubationStart = NaN;
                if ~isempty(AllHolds)
                    whichSegment = find(AllHolds>=0.9*thisTrialHolds(1),1,'first'); % any contiguous hold that would trigger offset perturbation
                    if isempty(whichSegment)
                        whichSegment = ...
                            find(cumsum(AllHolds)>=0.9*thisTrialHolds(2),1,'first'); % any aggregate hold that would trigger offset perturbation
                        if ~isempty(whichSegment)
                            deltaT = sum(AllHolds)-0.9*thisTrialHolds(2);
                            PertubationStart = SampleRate* ...
                                (TrialInfo.InZone{i}(whichSegment,2) - deltaT); % in indices w.r.t. TrialStart
                        end
                    else
                        PertubationStart = SampleRate* ...
                            (TrialInfo.InZone{i}(whichSegment,1) + thisTrialHolds(1)*0.9); % in indices w.r.t. TrialStart
                    end
                end
                PerturbationEvents(i,:) = [PertubationStart NaN NaN];
            else
                if ~isempty(TrialInfo.Perturbation{i,2})
                    PerturbationEvents(i,:) = TrialInfo.Perturbation{i,2};
                else
                    PerturbationEvents(i,:) = [NaN NaN NaN];
                end
            end
        end
        
    end
    
    if any(strcmp(TrialInfo.Perturbation(x,1),'RuleReversal'))
        load(MySession,'SampleRate');
        PerturbationEvents(1:nTrials,1) = 0;
        % mark all flipped trials as one
        PerturbationEvents(find(strcmp(TrialInfo.Perturbation(:,1),'RuleReversal')),1) = 1;
    end
end

if sniffwarpmethod
    % for aligning to sniffs
    % recompute Trialstart in behavior timebase - only need this because of
    % sample drops at trial start
    trialtimes = TrialInfo.SessionTimestamps(:,[1 2 1]);
    funkytrials = find(TrialInfo.TimeStampsDropped);
    if ~isempty(funkytrials)
        trialtimes(funkytrials,1) = trialtimes(funkytrials,2) - TTLs.Trial(funkytrials,3);
    end
    
    load(MySession,'SniffTS'); % sniff timestamps in behavior timebase
    
    if ~isempty(SniffTS)
        % get a cell array trial aligned sniffs
        for whichtrial = 1: nTrials % every trial
            if whichtrial==1
                first_inhalation = 1;
                last_inhalation  = find(SniffTS(:,1)<=trialtimes(whichtrial+1,1),1,'last');
            elseif whichtrial == nTrials
                first_inhalation = find(SniffTS(:,1)<=trialtimes(whichtrial-1,2),1,'last');
                last_inhalation  = size(SniffTS,1);
            else
                first_inhalation = find(SniffTS(:,1)<=trialtimes(whichtrial-1,2),1,'last');
                last_inhalation  = find(SniffTS(:,1)<=trialtimes(whichtrial+1,1),1,'last');
            end
            
            mysniffs = SniffTS(first_inhalation:last_inhalation,:) - trialtimes(whichtrial,1);
            mysniffs(:,4) = (1:size(mysniffs,1))' - numel(find(mysniffs(:,1)<=0));
            TrialAlignedSniffs{whichtrial} = mysniffs;
            
        end
    end
end

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
        
        %         if ~isempty(previousTrialspikes)
        %             previousTrialspikes(previousTrialspikes<-1) = [];
        %         end
        
        AlignedSpikes{whichtrial,whichunit} = {horzcat(previousTrialspikes',thisTrialspikes')};
        
        if sniffwarpmethod
            % convert every spike to sniff phase
            myspikes = horzcat(previousTrialspikes',thisTrialspikes');
            SniffAlignedSpikes{whichtrial,whichunit} = { WhichSniffPhase(myspikes,TrialAlignedSniffs{whichtrial},'warpmethod',sniffwarpmethod) };
        else
            SniffAlignedSpikes{whichtrial,whichunit} = { horzcat(previousTrialspikes',thisTrialspikes') };
        end
        
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
            
            if sniffwarpmethod
                % alignEventTimes to Sniffs
                Events_phase(whichtrial,:) = WhichSniffPhase(Events(whichtrial,:),TrialAlignedSniffs{whichtrial},'warpmethod',sniffwarpmethod);
                
                % also recalculate InZone periods in sniff base
                if ~isempty(TrialInfo.InZone{whichtrial})
                    TrialInfo.InZonePhase{whichtrial}(:,1) = WhichSniffPhase(TrialInfo.InZone{whichtrial}(:,1),TrialAlignedSniffs{whichtrial});
                    TrialInfo.InZonePhase{whichtrial}(:,2) = WhichSniffPhase(TrialInfo.InZone{whichtrial}(:,2),TrialAlignedSniffs{whichtrial});
                end
            else
                Events_phase(whichtrial,:) = Events(whichtrial,:);
            end
        end
    end
end

end