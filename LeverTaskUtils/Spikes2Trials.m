function [SingleUnits] = Spikes2Trials(SingleUnits, BehaviorTTLs, TuningTTLs)

%% function to label spike times by trials
% inputs: BehaviorTTLs - Trial On-Off times (offset corrected) in Oeps timebase

%% defaults
global startoffset; % = 1; % seconds

% for the first set of TTLs
for myUnit = 1:length(SingleUnits) % for each cluster
    
    allspikes = SingleUnits(myUnit).spikes; % in seconds
    SingleUnits(myUnit).trialtags = NaN*allspikes;
    SingleUnits(myUnit).trialalignedspikes = allspikes;
    
    % Assign spikes to trials
    % Spikes preceding trial one get Id = 0
    % all other Spikes belong to trial that started just before
    
    % behavior TTLs
    for mytrial = 1:size(BehaviorTTLs,1)
        tstart = BehaviorTTLs(mytrial,1);
        if mytrial == 1
            SingleUnits(myUnit).trialtags(find(allspikes<tstart)) = 0;
        end
        if mytrial < size(BehaviorTTLs,1)
            tstop = BehaviorTTLs(mytrial+1,1);
        else
            tstop = BehaviorTTLs(mytrial,2) + startoffset;
        end
        SingleUnits(myUnit).trialtags(find(allspikes>=tstart & allspikes<tstop)) = mytrial;
        SingleUnits(myUnit).trialalignedspikes(find(allspikes>=tstart & allspikes<tstop)) = ...
            SingleUnits(myUnit).spikes(find(allspikes>=tstart & allspikes<tstop)) - tstart;
    end
    
    % Tuning TTLs - including passive replays
    if nargin>2
        if ~isempty(TuningTTLs)
            for mytrial = 1:size(TuningTTLs,1)
                tstart = TuningTTLs(mytrial,1);
                if mytrial < size(TuningTTLs,1)
                    tstop = TuningTTLs(mytrial+1,1);
                else
                    tstop = TuningTTLs(mytrial,2) + startoffset;
                end
                SingleUnits(myUnit).trialtags(find(allspikes>=tstart & allspikes<tstop)) = -TuningTTLs(mytrial,end);
                SingleUnits(myUnit).trialalignedspikes(find(allspikes>=tstart & allspikes<tstop)) = ...
                    SingleUnits(myUnit).spikes(find(allspikes>=tstart & allspikes<tstop)) - tstart;
            end
        end
    end
end

end
