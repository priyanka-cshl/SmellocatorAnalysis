function [AlignedSpikes, TetrodeOrder] = TrialAlignedSpikeTimes_Tuning(SingleUnits,TTLs)

N = size(SingleUnits,2); % total units
nTrials = size(TTLs,1);

for whichunit = 1:N % every unit
    thisUnitspikes = SingleUnits(whichunit).trialalignedspikes;
    trialtags = SingleUnits(whichunit).trialtags;
    TetrodeOrder(whichunit,:) = [SingleUnits(whichunit).tetrode SingleUnits(whichunit).id]; % tetrode and phy cluster ID
    for whichtrial = 1: nTrials % every trial 
        original_trial = TTLs(whichtrial,8); % actual trial # (including behavior trials)
        thisTrialspikes = thisUnitspikes(trialtags==-original_trial);
        if whichtrial>1
            previousTrialspikes = thisUnitspikes(trialtags==-(original_trial-1)) - ...
                (TTLs(whichtrial,1) - TTLs(whichtrial-1,1));
            
            % ignore spikes that happened before previous trial off
            prevToff = TTLs(whichtrial-1,2) - TTLs(whichtrial,1);
            if ~isempty(previousTrialspikes)
                previousTrialspikes(previousTrialspikes<prevToff) = [];
            end
        else
            previousTrialspikes = thisUnitspikes(trialtags==-(original_trial-1)) - ...
                (TTLs(whichtrial,1) - 0);
        end
        
%         if ~isempty(previousTrialspikes)
%             previousTrialspikes(previousTrialspikes<-1) = [];
%         end

        AlignedSpikes{whichtrial,whichunit} = {horzcat(previousTrialspikes',thisTrialspikes')};
    end
end

end