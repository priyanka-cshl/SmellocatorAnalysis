function [TrialAlignedSniffs, SniffAlignedSpikes, ReplayInfo] = ...
    SniffAlignedSpikeTimes_Replays(SingleUnits,TTLs,ReplayTTLs,TrialInfo,OpenLoop,MySession,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('sniffwarpmethod', 3, @(x) isnumeric(x)); % 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency

% extract values from the inputParser
params.parse(varargin{:});
sniffwarpmethod = params.Results.sniffwarpmethod;

N = size(SingleUnits,2); % total units

% there are multiple templates, each template should have a perturbed trial
% for each template, pull out the odor sequence and identify which trial was perturbed
nTemplates = size(OpenLoop.TemplateTraces.TrialIDs,2);
Templates.Trials = [];
Templates.Odors = [];
for i = 1:nTemplates
    whichtrials = OpenLoop.TemplateTraces.TrialIDs{i};
    Templates.Trials(i,1:numel(whichtrials)) = whichtrials;
    Templates.Odors(i,1:numel(whichtrials))  = TrialInfo.Odor(whichtrials);
    f = find(~strcmp(TrialInfo.Perturbation(whichtrials,1),'OL-Template'));
    if ~isempty(f)
        Templates.Trials(i,f(1)) = -Templates.Trials(i,f(1));
    end
end

nTrials = numel(ReplayTTLs.TrialID); % how many replays
subtrials = [];
for j = 1:nTrials
    subtrials = vertcat(subtrials,size(ReplayTTLs.OdorValve{j},1));
end
nsubtrials = mode(subtrials);

% load sniffs
load(MySession,'SniffTS_passive', 'Passive_Timestamp_adjust'); % sniff timestamps in behavior timebase
SniffTS_passive(:,1:3) = SniffTS_passive(:,1:3) + Passive_Timestamp_adjust; % Sniff Times in OEPS timebase

if ~isempty(SniffTS_passive)
    % get trial times for all replay trials?
    trialcount = 0;
    for whichReplay = 1: nTrials % every replay
        thisReplay = ReplayTTLs.TrialID(whichReplay); % Replay Trial ID
        OdorTTLs = ReplayTTLs.OdorValve{whichReplay}; % TS of odor valve ON-OFF w.r.t. Trial start TTL
        
        if size(OdorTTLs,1) >= nsubtrials % valid replay
            % which template does this replay correspond to
            [~, whichtemplate] = ismember(OdorTTLs(:,4)',Templates.Odors(:,1:numel(OdorTTLs(:,4)')), 'rows');
            
            templatetrials = abs(Templates.Trials(whichtemplate(1),:));
            FirstTrialDuration = TrialInfo.Duration(templatetrials(1));
            perturbedtrial = find(Templates.Trials(whichtemplate(1),:)<0);
            
            if isempty(perturbedtrial)
                perturbedtrial = NaN;
            end
            
            % Add one more column for TrialStart w.r.t. Odor ON
            OdorTTLs(:,end+1) = OdorTTLs(:,1) - TrialInfo.OdorStart(templatetrials(find(templatetrials)),1);
            
            % force first subtrial to start at ~0
            OdorTTLs(1,end) = OdorTTLs(1,2) - FirstTrialDuration;
            
            % convert TS to real OEPS time base
            OdorTTLs(:,[1 2 5]) = OdorTTLs(:,[1 2 5]) + TTLs.Trial(thisReplay,1);
            
            for j = 1:nsubtrials % every trial within a replay stretch
                trialcount = trialcount + 1;
                
                if j == 1
                    t1 = TTLs.Trial(thisReplay-1,2);
                else
                    t1 = OdorTTLs(j-1,2); % previous trial's odor OFF (OEPS) = Trial OFF
                end
                
                t2 = OdorTTLs(j,5); % Trial ON (OEPS)
                t3 = OdorTTLs(j,2); % Trial OFF (OEPS) - same as odor OFF
                
                % t1 and t2 are previous trial OFF, this trial start and stop times in OEPS timebase
                % sniffs are already in OEPS timebase
                first_inhalation = find(SniffTS_passive(:,1)>=t1,1,'first'); % after prev trial OFF
                last_inhalation  = find(SniffTS_passive(:,2)<=t3,1,'last'); % before trial end
                
                mysniffs = SniffTS_passive(first_inhalation:last_inhalation,1:3) - t2; % -ve timestamps are in ITI
                mysniffs(:,4) = (1:size(mysniffs,1))' - numel(find(mysniffs(:,1)<=0));
                mysniffs(:,6) = SniffTS_passive(first_inhalation:last_inhalation,4); % odor location
                
                % flag sniffs where entire inhalation period was in the target
                % zone?
                whichtrial = templatetrials(j);
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
                    SniffTS(first_inhalation-1:last_inhalation-1,1:3) - t2, ...
                    SniffTS(first_inhalation-1:last_inhalation-1,4), ...
                    SniffTS(first_inhalation+1:last_inhalation+1,1:3) - t2, ...
                    SniffTS(first_inhalation+1:last_inhalation+1,4) );
                
                TrialAlignedSniffs{trialcount} = mysniffs;
                
                % get trial aligned spikes
                for whichunit = 1:N % every unit
                    thisUnitspikes  = SingleUnits(whichunit).spikes;
                    thisTrialspikes = thisUnitspikes(find(thisUnitspikes>=t1,1,'first'):find(thisUnitspikes<=t3,1,'last'));
                    thisTrialspikes = thisTrialspikes - t2; % spiketimes w.r.t. trial start
                    
                    % convert every spike to sniff phase
                    SniffAlignedSpikes{trialcount,whichunit} = { WhichSniffPhase(thisTrialspikes,TrialAlignedSniffs{trialcount},'warpmethod',sniffwarpmethod) };
                    
                end
                
                % Useful stuff for parsing trials into odors and perturbations
                ReplayInfo.Odor(trialcount,1)           = TrialInfo.Odor(templatetrials(j));
                ReplayInfo.Perturbed(trialcount,1)      = isequal(j,perturbedtrial);
                ReplayInfo.Perturbation(trialcount,:)   = TrialInfo.Perturbation(templatetrials(j),:);
                if thisReplay <= TrialInfo.TrialID(end)
                    ReplayInfo.TrialID(trialcount)          =  whichReplay + j/100;
                else
                    ReplayInfo.TrialID(trialcount)          = -(whichReplay + j/100);
                end
                
            end
        end
    end
end