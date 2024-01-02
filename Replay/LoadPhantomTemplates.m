function [Templates, TrialInfo] = LoadPhantomTemplates(MyPhantomSession)
%% a hack code for loading the templates for passive replay from another mouse
% such that replay trials can still be made sense of - in terms of chunking
% trials etc

[TracesOut, ColNames, TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MyPhantomSession); % LoadProcessedSession; % loads relevant variables

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

end