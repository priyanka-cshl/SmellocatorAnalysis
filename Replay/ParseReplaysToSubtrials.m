function [ReplayInfo] = ParseReplaysToSubtrials(OpenLoop,TrialInfo,ReplayTTLs,TTLs)

%   function to parse out the passive replays into subtrials
%   and get subtrial params

%% STEP 1: Is there one (OpenLoop Sessions) 
% or many (Passive perturbation replays) templates?
nTemplates = size(OpenLoop.TemplateTraces.TrialIDs,2);

% for each template, pull out the odor sequence, and    
% if nTemplates > 1
%   perturbation containing passive replays
%   each template should have a perturbed trial
%   and identify which trial was perturbed
    
Templates.Trials = [];
Templates.Odors = [];

for i = 1:nTemplates
    % every template is composed of many subtrials
    % find the trial indices of those subtrials
    whichtrials = OpenLoop.TemplateTraces.TrialIDs{i};
    Templates.Trials(i,1:numel(whichtrials)) = whichtrials; % indices
    Templates.Odors(i,1:numel(whichtrials))  = TrialInfo.Odor(whichtrials); % odor identities
    
    % perturbed subtrials are flagged with the name of the perturbation
    % other falnking subtrials are tagged as 'OL-Template'
    % find the perturbed trial and make its trial index negative
    f = find(~strcmp(TrialInfo.Perturbation(whichtrials,1),'OL-Template'));
    if ~isempty(f)
        Templates.Trials(i,f(1)) = -Templates.Trials(i,f(1));
    end
end

%% STEP 2: Parse the individual passive relays (one long trial) into subtrials
% using info from the template trials
trialcount  = 0;
ReplayInfo  = [];
nReplays    = numel(ReplayTTLs.TrialID); % how many replays
if nTemplates > 1
    nSubtrials  = mode( cellfun(@(x) size(x,1), ReplayTTLs.OdorValve)); % typical no. of trials in a template
else
    nSubtrials  = min( cellfun(@(x) size(x,1), ReplayTTLs.OdorValve));
end
    
for whichReplay = 1:nReplays % every replay
    thisReplay = ReplayTTLs.TrialID(whichReplay); % Replay Trial ID
    OdorTTLs = ReplayTTLs.OdorValve{whichReplay}; % TimeStamp of odor valve ON-OFF w.r.t. Trial start TTL - Behavior base and odor identity
    
    if size(OdorTTLs,1) >= nSubtrials % valid replay
        % get the odor sequence for this Replay's subtrials from OdorTTLs
        % sometimes there is an extra odor trial in the beginning 
        % because of quick transiion from prev trial's odor to the new one
        % - delete that
        if OdorTTLs(1,2) < 0.1 
            OdorTTLs(1,:) = [];
        end
        thisReplayOdors = OdorTTLs(:,4)'; % [1 x nSubtrials]

        % find the matching template for this replay using the odor sequence
        try
            [~, whichtemplate] = ismember(thisReplayOdors,Templates.Odors(:,1:numel(thisReplayOdors)), 'rows');
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
        
        % flag the perturbed trial
        perturbedtrial = find(Templates.Trials(whichtemplate(1),:)<0);
        if isempty(perturbedtrial)
            perturbedtrial = NaN;
        end
        
        % ------------------------------------------------------
        % After Matching template is found
        % ------------------------------------------------------
        %% STEP 2.1: compute SubTrial Start/Stop times for the Replay using OdorTTLs
        % Trial Start = Odor Start - Pre-Trial OdorON duration
        % Trial Stop = Odor Stop
        
        % get subtrial indices from the matching Template
        templatetrials = abs(Templates.Trials(whichtemplate(1),:)); 
        
        % Add a column in OdorTTLs (Col 5) for TrialStart w.r.t. Odor ON
        OdorTTLs(:,5) = OdorTTLs(:,1) - TrialInfo.OdorStart(templatetrials,1);
        % OdorTTLs(:,end+1) = OdorTTLs(:,1) - TrialInfo.OdorStart(templatetrials(find(templatetrials)),1);
        
        % get duration of the first trial - for aligning the first trial
        % (again due to quick odor transitions in the beginning)
        FirstTrialDuration = TrialInfo.Duration(templatetrials(1));
        % force first subtrial to start at ~0
        OdorTTLs(1,5) = OdorTTLs(1,2) - FirstTrialDuration;
        
        % convert computed TS from behavior to real OEPS time base
        OdorTTLs(:,[1 2 5]) = OdorTTLs(:,[1 2 5]) + TTLs.Trial(thisReplay,1);
        % Col 3 is trial duraion and 4 is odor identity
        
        %% STEP 2.2: Compile a list of Trial On/Off times for the current subtrial and also the previous trial's OFF and next Trial's Start
        for j = 1:nSubtrials % every trial within a replay stretch
            trialcount = trialcount + 1;
            ReplayInfo.trialtimes(trialcount,1:2) = OdorTTLs(j,[5 2]); % Trial ON and OFF (OEPS base)
            % previous Trial's OFF
            if j == 1
                ReplayInfo.trialtimes(trialcount,3) = TTLs.Trial(thisReplay-1,2);
            else
                ReplayInfo.trialtimes(trialcount,3) = OdorTTLs(j-1,2); % previous trial's odor OFF (OEPS) = Trial OFF
            end
            % next Trial's ON
            if j<nSubtrials
                ReplayInfo.trialtimes(trialcount,4) = OdorTTLs(j+1,1); % next odor ON (OEPS)
            else
                ReplayInfo.trialtimes(trialcount,4) = TTLs.Trial(thisReplay+1,1);
            end
            
            % Useful stuff for parsing trials into odors and perturbations
            ReplayInfo.Odor(trialcount,1)           = TrialInfo.Odor(templatetrials(j));
            ReplayInfo.OdorStart(trialcount,1)      = TrialInfo.OdorStart(templatetrials(j));
            ReplayInfo.Perturbed(trialcount,1)      = isequal(j,perturbedtrial);
            ReplayInfo.Perturbation(trialcount,:)   = TrialInfo.Perturbation(templatetrials(j),:);
            if thisReplay <= TrialInfo.TrialID(end)
                ReplayInfo.TrialID(trialcount,1)    =  whichReplay + j/100;
            else
                ReplayInfo.TrialID(trialcount,1)    = -(whichReplay + j/100);
            end
            ReplayInfo.OriginalTrialID(trialcount,1)= thisReplay; 
            ReplayInfo.TemplateInfo(trialcount,1:2)   = [whichtemplate templatetrials(j)];
        end
        
    end
end

end