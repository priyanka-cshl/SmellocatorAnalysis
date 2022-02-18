function [AlignedSpikes, TuningTTLs] = TrialAlignedSpikeTimesPassiveTuning(SingleUnits,TTLs,TuningTTLs)

N = size(SingleUnits,2); % total units

% remove any passive replay trials
TuningTTLs(find(isnan(TuningTTLs(:,5))),:) = [];
% remove any trials for which motor location canot be confirmed
TuningTTLs(find(isnan(TuningTTLs(:,7))),:) = [];

for whichunit = 1:N % every unit
    thisUnitspikes = SingleUnits(whichunit).trialalignedspikes;
    trialtags = SingleUnits(whichunit).trialtags;
    
    for thistrial = 1: size(TuningTTLs,1) % every trial 
        whichtrial = TuningTTLs(thistrial,end); % trial tag as per TTLs
        thisTrialspikes = thisUnitspikes(abs(trialtags)==whichtrial);
        if whichtrial>1
            previousTrialspikes = thisUnitspikes(abs(trialtags)==whichtrial-1) - ...
                (TTLs.Trial(whichtrial,1) - TTLs.Trial(whichtrial-1,1));
            
            % ignore spikes that happened before previous trial off
            prevToff = TTLs.Trial(whichtrial-1,2) - TTLs.Trial(whichtrial,1);
            if ~isempty(previousTrialspikes)
                previousTrialspikes(previousTrialspikes<prevToff) = [];
            end
        else
            previousTrialspikes = thisUnitspikes(abs(trialtags)==whichtrial-1) - ...
                (TTLs.Trial(whichtrial,1) - 0);
        end
        
        AlignedSpikes{thistrial,whichunit} = {horzcat(previousTrialspikes',thisTrialspikes')};
    end
end

end