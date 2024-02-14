function [EventAlignedSpikes] = EventAlignedSpikeTimes(SingleUnits,EventWindows)
% function to get spikes aligned to a particular event window eg. trial
% Input: SingleUnits (data struct with unitwise spike times in OEPS base)
%        EventWindows [n x 1-4] in OEPS time base
%           n = no. of events
%           Cols: [WindowStart WindowEnd PrevWindowEnd NextWindowStart] 
%           If nCols = 1, use default window sizes
% Output : EventAlignedSpikes{nEvents,nUnits}
%           spiketimes are calculated w.r.t. WindowStart

defaultWindow = 1; % seconds
if size(EventWindows,2) < 4
    EventWindows(:,3) = EventWindows(:,1) - defaultWindow;
    if size(EventWindows,2) < 2
        EventWindows(:,[2 4]) = EventWindows(:,1) + defaultWindow;
    elseif size(EventWindows,2) < 3
        EventWindows(:,4) = EventWindows(:,2) + defaultWindow;
    end
end

%% Align Spikes to EventWindows
for whichunit = 1:size(SingleUnits,2) % every unit
    thisUnitspikes = SingleUnits(whichunit).spikes; % raw spiketimes in OEPS unit

    for whichEvent = 1: size(EventWindows,1) % every Event (trial)
        thisTrialspikes = thisUnitspikes((thisUnitspikes>EventWindows(whichEvent,3))&&(thisUnitspikes<=EventWindows(whichEvent,4)));
        thisTrialspikes = thisTrialspikes - EventWindows(whichEvent,1);
        EventAlignedSpikes{whichEvent,whichunit} = { thisTrialspikes };
        
        % convert every spike to sniff phase
        SniffAlignedSpikes{whichEvent,whichunit} = { WhichSniffPhase(myspikes,TrialAlignedSniffs{whichEvent},'warpmethod',sniffwarpmethod) };
        
        if whichunit == 1
            x1 = TTLs.Trial(whichEvent,1);
            x2 = TTLs.Trial(whichEvent,2) + 0.2; % next trial cannot start sooner than that
            reward = TTLs.Reward(find((TTLs.Reward(:,1)>x1)&(TTLs.Reward(:,1)<x2)),1);
            if isempty(reward)
                reward = NaN;
            end
            % OdorStart, Reward, TrialOFF, First TZ entry
            Events(whichEvent,1:4) = [TTLs.Trial(whichEvent,4), ...
                (reward(1) - x1), ...
                TTLs.Trial(whichEvent,3), ...
                TrialInfo.TargetEntry(whichEvent,1)];
            
            if ~isempty(PerturbationEvents)
                Events(whichEvent,5) = PerturbationEvents(whichEvent,1)/SampleRate;
            end
            
            % alignEventTimes to Sniffs
            Events_phase(whichEvent,:) = WhichSniffPhase(Events(whichEvent,:),TrialAlignedSniffs{whichEvent},'warpmethod',sniffwarpmethod);
            
            % also recalculate InZone periods in sniff base
            if ~isempty(TrialInfo.InZone{whichEvent})
                TrialInfo.InZonePhase{whichEvent}(:,1) = WhichSniffPhase(TrialInfo.InZone{whichEvent}(:,1),TrialAlignedSniffs{whichEvent});
                TrialInfo.InZonePhase{whichEvent}(:,2) = WhichSniffPhase(TrialInfo.InZone{whichEvent}(:,2),TrialAlignedSniffs{whichEvent});
            end

        end
        
    end
end

end