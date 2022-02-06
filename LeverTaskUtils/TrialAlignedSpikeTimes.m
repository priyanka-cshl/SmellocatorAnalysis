function [AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,nTrials,TrialInfo)

N = size(SingleUnits,2); % total units
for whichunit = 1:N % every unit
    thisUnitspikes = SingleUnits(whichunit).trialalignedspikes;
    trialtags = SingleUnits(whichunit).trialtags;
    
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
            Events(whichtrial,:) = [TTLs.Trial(whichtrial,4), ...
                                    (reward(1) - x1), ...
                                    TTLs.Trial(whichtrial,3), ...
                                    TrialInfo.TargetEntry(whichtrial,1)];
        end
    end
end

end