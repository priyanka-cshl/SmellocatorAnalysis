function [ReplayAlignedSpikes, ReplayEvents, ReplayInfo] = PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,ReplayTTLs,TrialInfo,Events,OpenLoop)

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

N = size(SingleUnits,2); % total units
for whichunit = 1:N % every unit
    %thisUnitspikes = SingleUnits(whichunit).trialalignedspikes;
    thisUnitspikes = SingleUnits(whichunit).spikes;
    %trialtags = abs(SingleUnits(whichunit).trialtags); % tuning TTLs are in negative
    trialcount = 0;
    for whichReplay = 1: nTrials % every replay
        %templatetrials = find(strcmp(TrialInfo.Perturbation,'OL-Template'));
        %FirstTrialDuration = TrialInfo.Duration(templatetrials(1));
        
        thisReplay = ReplayTTLs.TrialID(whichReplay); % Replay Trial ID
        OdorTTLs = ReplayTTLs.OdorValve{whichReplay}; % TS of odor valve ON-OFF w.r.t. Trial start TTL
        
        if size(OdorTTLs,1) >= nsubtrials
%             if size(OdorTTLs,1)>nsubtrials % some replays have an extra odor transition at the very beginning - to shut off odor from pre-replay trial
%                 OdorTTLs(1,:) = [];
%             end
            
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
                
                if j<nsubtrials'
                    t3 = OdorTTLs(j+1,1); % next odor ON (OEPS)
                else
                    t3 = TTLs.Trial(thisReplay+1,1);
                end
                
                thisReplaySpikes = thisUnitspikes((thisUnitspikes>t1)&(thisUnitspikes<t3)) - t2;
                ReplayAlignedSpikes{trialcount,whichunit} = {thisReplaySpikes};
                
                if whichunit == 1
                    
                    ReplayEvents(trialcount,2:5)            = Events(templatetrials(j),2:5);
                    ReplayEvents(trialcount,1)              = OdorTTLs(j,1) - OdorTTLs(j,5); % odor ON w.r.t. trial start
                    ReplayInfo.Odor(trialcount,1)           = TrialInfo.Odor(templatetrials(j));
                    ReplayInfo.TargetZoneType(trialcount,1) = TrialInfo.TargetZoneType(templatetrials(j));
                    ReplayInfo.Duration(trialcount,1)       = TrialInfo.Duration(templatetrials(j));
                    ReplayInfo.InZone{trialcount}           = TrialInfo.InZone{templatetrials(j)};
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
end

end