function [TrialAlignedSniffs, SniffAlignedSpikes, TrialAlignedSpikes, TetrodeOrder] = SniffAlignedSpikeTimes(SingleUnits,TTLs,nTrials,TrialInfo,MySession,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('sniffwarpmethod', 3, @(x) isnumeric(x)); % 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency

% extract values from the inputParser
params.parse(varargin{:});
sniffwarpmethod = params.Results.sniffwarpmethod;

N = size(SingleUnits,2); % total units

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
            first_inhalation = 2;
            last_inhalation  = find(SniffTS(:,2)<=trialtimes(whichtrial,2),1,'last');
        else
            first_inhalation = find(SniffTS(:,1)>=trialtimes(whichtrial-1,2),1,'first'); % after prev trial OFF
            last_inhalation  = find(SniffTS(:,2)<=trialtimes(whichtrial,2),1,'last'); % before trial end
        end
        
        mysniffs = SniffTS(first_inhalation:last_inhalation,1:3) - trialtimes(whichtrial,1); % -ve timestamps are in ITI
        mysniffs(:,4) = (1:size(mysniffs,1))' - numel(find(mysniffs(:,1)<=0));
        mysniffs(:,6) = SniffTS(first_inhalation:last_inhalation,4); % odor location 
        
        % flag sniffs where entire inhalation period was in the target
        % zone?
        if ~isempty(TrialInfo.InZone{whichtrial})
            
            whichsniffs = find(mysniffs(:,1)>=TrialInfo.TargetEntry(whichtrial,1));
            mysniffs(whichsniffs,5) = 1;
            for entries = 1:size(TrialInfo.InZone{whichtrial},1)
                stay = TrialInfo.InZone{whichtrial}(entries,:);
                whichsniffs = intersect( (find(mysniffs(:,1)>=stay(1))) , ...
                    (find(mysniffs(:,2)<=stay(2))) );
                mysniffs(whichsniffs,5) = 2;
            end
        end
        mysniffs(mysniffs(:,3)<=0,5) = -1;
        
        % also note the sniff before and sniff after
        mysniffs = horzcat(mysniffs, ...
            SniffTS(first_inhalation-1:last_inhalation-1,1:3) - trialtimes(whichtrial,1), ...
            SniffTS(first_inhalation-1:last_inhalation-1,4), ...
            SniffTS(first_inhalation+1:last_inhalation+1,1:3) - trialtimes(whichtrial,1), ...
            SniffTS(first_inhalation+1:last_inhalation+1,4) );
        
        TrialAlignedSniffs{whichtrial} = mysniffs;
        
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
                
        % convert every spike to sniff phase
        myspikes = horzcat(previousTrialspikes',thisTrialspikes');
        TrialAlignedSpikes{whichtrial,whichunit} = { myspikes };
        SniffAlignedSpikes{whichtrial,whichunit} = { WhichSniffPhase(myspikes,TrialAlignedSniffs{whichtrial},'warpmethod',sniffwarpmethod) };
    end
end

end