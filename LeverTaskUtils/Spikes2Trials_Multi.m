function [SingleUnits, TrialStretch] = Spikes2Trials_Multi(SingleUnits, allTTLs, BehaviorTrials, TuningTTLs_all)

%% function to label spike times by trials
% inputs: TTLs - Trial On-Off times (offset corrected) in Oeps timebase

BehaviorTTLs = [];
for n = 1:size(BehaviorTrials,1) % no. of closed loop sessions
    BehaviorTTLs = [BehaviorTTLs; ...
                    horzcat(allTTLs(BehaviorTrials(n,1):BehaviorTrials(n,2),:), ...
                                    (BehaviorTrials(n,1):BehaviorTrials(n,2))') ];
                           %  (n-1)*1000+(BehaviorTrials(n,1):BehaviorTrials(n,2))') ];
end

%% defaults
global startoffset; % = 1; % seconds
TuningTrials = [];

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
        %SingleUnits(myUnit).trialtags(find(allspikes>=tstart & allspikes<tstop)) = mytrial;
        SingleUnits(myUnit).trialtags(find(allspikes>=tstart & allspikes<tstop)) = BehaviorTTLs(mytrial,end);
        SingleUnits(myUnit).trialalignedspikes(find(allspikes>=tstart & allspikes<tstop)) = ...
            SingleUnits(myUnit).spikes(find(allspikes>=tstart & allspikes<tstop)) - tstart;
    end
    
    % Tuning TTLs - including passive replays
    if nargin>2
        for n = 1:size(TuningTTLs_all,2) % no. of passive tuning sessions
            TuningTTLs = TuningTTLs_all{n};
            if ~isempty(TuningTTLs)
                for mytrial = 1:size(TuningTTLs,1)
                    tstart = TuningTTLs(mytrial,1);
                    if mytrial < size(TuningTTLs,1)
                        tstop = TuningTTLs(mytrial+1,1);
                    else
                        tstop = TuningTTLs(mytrial,2) + startoffset;
                    end
                    SingleUnits(myUnit).trialtags(find(allspikes>=tstart & allspikes<tstop)) = -TuningTTLs(mytrial,8); % original Trial ID
                    SingleUnits(myUnit).trialalignedspikes(find(allspikes>=tstart & allspikes<tstop)) = ...
                        SingleUnits(myUnit).spikes(find(allspikes>=tstart & allspikes<tstop)) - tstart;
                end
                if myUnit == 1
                    TuningTrials(n,2) = TuningTTLs(1,8);
                    TuningTrials(n,3) = TuningTTLs(mytrial,8);
                    TuningTrials(n,1) = -1;
                end
            end
        end
    end
end

% Send a concatenated trial idx list out for analysis later
TrialStretch = [ones(size(BehaviorTrials,1),1) BehaviorTrials];
TrialStretch = [TrialStretch; TuningTrials];

end
