function [ReplayAlignedSpikes, ReplayEvents, ReplayInfo, SniffAlignedReplaySpikes, ReplayEvents_phase] = ...
    PerturbationReplayAlignedSpikeTimes_v2(SingleUnits,TTLs,ReplayTTLs,TrialInfo,Events,OpenLoop,MySession,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('sniffwarpmethod', 0, @(x) isnumeric(x)); % 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency

% extract values from the inputParser
params.parse(varargin{:});
sniffwarpmethod = params.Results.sniffwarpmethod;

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

% % hack for S3_20230327_r0_processed.mat 
% % ---- managed to do halt replays from S1 instead 
% if contains(MySession,'S3_20230327_r0_processed.mat')
%     MyPhantomSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S1/S1_20230327_r0_processed.mat';
%     [Templates, TrialInfo] = LoadPhantomTemplates(MyPhantomSession);
% end

nTrials = numel(ReplayTTLs.TrialID); % how many replays
subtrials = [];
for j = 1:nTrials
    subtrials = vertcat(subtrials,size(ReplayTTLs.OdorValve{j},1));
end
nsubtrials = mode(subtrials);

if sniffwarpmethod
    load(MySession,'SniffTS_passive', 'Passive_Timestamp_adjust'); % sniff timestamps in behavior timebase
    SniffTS_passive(:,1:3) = SniffTS_passive(:,1:3) + Passive_Timestamp_adjust; % Sniff Times in OEPS timebase
end

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
            try
                [~, whichtemplate] = ismember(OdorTTLs(:,4)',Templates.Odors(:,1:numel(OdorTTLs(:,4)')), 'rows');
            catch
                disp(['no template match for replay# ',num2str(whichReplay)]);
                keyboard;
                continue;
            end
            
            if ~whichtemplate
                disp(['no template match for replay# ',num2str(whichReplay)]);
                keyboard;
                continue;
            end
                        
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
                
                if j<nsubtrials
                    t3 = OdorTTLs(j+1,1); % next odor ON (OEPS)
                else
                    t3 = TTLs.Trial(thisReplay+1,1);
                end
                
                thisReplaySpikes = thisUnitspikes((thisUnitspikes>t1)&(thisUnitspikes<t3)) - t2;
                ReplayAlignedSpikes{trialcount,whichunit} = {thisReplaySpikes};
                
                if sniffwarpmethod
                    thisReplaySniifs = SniffTS_passive(((SniffTS_passive(:,1)>=t1)&(SniffTS_passive(:,1)<=t3)),:) - t2;
                    thisReplaySniifs(:,4) = (1:size(thisReplaySniifs,1))' - numel(find(thisReplaySniifs(:,1)<=0));
                    
                    % convert every spike to sniff phase
                    SniffAlignedReplaySpikes{trialcount,whichunit} = { WhichSniffPhase(thisReplaySpikes,thisReplaySniifs,'warpmethod',sniffwarpmethod) };
                else
                    SniffAlignedReplaySpikes{trialcount,whichunit} = { thisReplaySpikes };
                end
                
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
                    
                    if sniffwarpmethod
                        % alignEventTimes to Sniffs
                        ReplayEvents_phase(trialcount,:) = WhichSniffPhase(ReplayEvents(trialcount,:),thisReplaySniifs,'warpmethod',sniffwarpmethod);
                        
                        % also recalculate InZone periods in sniff base
                        if ~isempty(ReplayInfo.InZone{trialcount})
                            ReplayInfo.InZonePhase{trialcount}(:,1) = WhichSniffPhase(ReplayInfo.InZone{trialcount}(:,1),thisReplaySniifs);
                            ReplayInfo.InZonePhase{trialcount}(:,2) = WhichSniffPhase(ReplayInfo.InZone{trialcount}(:,2),thisReplaySniifs);
                        end
                    else
                        ReplayEvents_phase(trialcount,:) = Events(trialcount,:);
                    end
                end
            end
        end
    end
end

end