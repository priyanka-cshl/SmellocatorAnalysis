function [AlignedSpikes, Events, TetrodeOrder] = TrialAlignedSpikeTimes(SingleUnits,TTLs,nTrials,TrialInfo,MySession)

N = size(SingleUnits,2); % total units

PerturbationEvents = [];
if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    
    if any(strcmp(TrialInfo.Perturbation(x,1),'Halt-Flip')) || ...
            any(strcmp(TrialInfo.Perturbation(x,1),'Halt-Flip-Template'))
        
        HaltTrials = [find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')); ...
                      find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip-Template')) ];
                  
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
                PerturbationEvents(i,:) = TrialInfo.Perturbation{i,2};
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

for whichunit = 1:N % every unit
    thisUnitspikes = SingleUnits(whichunit).trialalignedspikes;
    trialtags = SingleUnits(whichunit).trialtags;
    TetrodeOrder(whichunit) = SingleUnits(whichunit).tetrode;
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
        end
    end
end

end