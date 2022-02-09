function [ReplayAlignedSpikes, ReplayEvents, ReplayInfo] = ReplayAlignedSpikeTimes(SingleUnits,TTLs,ReplayTTLs,TrialInfo,Events)

nTrials = numel(ReplayTTLs.TrialID);
subtrials = [];
for j = 1:nTrials
    subtrials = vertcat(subtrials,size(ReplayTTLs.OdorValve{j},1));
end
nsubtrials = min(subtrials);

templatetrials = find(strcmp(TrialInfo.Perturbation,'OL-Template'));

N = size(SingleUnits,2); % total units
for whichunit = 1:N % every unit
    thisUnitspikes = SingleUnits(whichunit).trialalignedspikes;
    trialtags = abs(SingleUnits(whichunit).trialtags); % tuning TTLs are in negative
    trialcount = 0;
    for whichtrial = 1: nTrials % every trial 
        thisTrial = ReplayTTLs.TrialID(whichtrial);
        thisTrialspikes = thisUnitspikes(trialtags==thisTrial);
        previousTrialspikes = thisUnitspikes(trialtags==thisTrial-1) - ...
            (TTLs.Trial(thisTrial,1) - TTLs.Trial(thisTrial-1,1));
        thisTrialspikes = horzcat(previousTrialspikes',thisTrialspikes');
        
        subTTLs = ReplayTTLs.OdorValve{whichtrial};
        if size(subTTLs,1)>nsubtrials
            subTTLs(1,:) = [];
        end
            
        for j = 1:nsubtrials
            trialcount = trialcount + 1;
            
            if j == 1
                t1 = min(thisTrialspikes)-1; 
                t2 = 0; % Trial Start
            else
                t1 = subTTLs(j-1,2); % previous trial's odor OFF (OEPS) = Trial OFF
                t2 = subTTLs(j,1); % odor ON (OEPS)
                t2 = t2 - TrialInfo.OdorStart(templatetrials(j),1); % convert to trial start 
            end
            
            if j<nsubtrials
                t3 = subTTLs(j+1,1); % next trial ON (OEPS)
            else
                t3 = TTLs.Trial(thisTrial+1,1);
            end
            
            thisReplaySpikes = thisTrialspikes((thisTrialspikes>t1)&(thisTrialspikes<t3)) - t2;
            ReplayAlignedSpikes{trialcount,whichunit} = {thisReplaySpikes};
            
            if whichunit == 1
                ReplayEvents(trialcount,1:4) = Events(templatetrials(j),:);
                ReplayInfo.Odor(trialcount,1)           = TrialInfo.Odor(templatetrials(j));
                ReplayInfo.TargetZoneType(trialcount,1) = TrialInfo.TargetZoneType(templatetrials(j));
                ReplayInfo.Duration(trialcount,1)       = TrialInfo.Duration(templatetrials(j));
                ReplayInfo.InZone{trialcount}           = TrialInfo.InZone{templatetrials(j)};
                if thisTrial <= TrialInfo.TrialID(end)
                    ReplayInfo.TrialID(trialcount)          =  whichtrial + j/100;
                else
                    ReplayInfo.TrialID(trialcount)          = -(whichtrial + j/100);
                end
            end
        end
    end
end

end