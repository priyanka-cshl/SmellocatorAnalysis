function [] = ProcessReplays(Traces, TrialInfo, TTLs, ReplayTTLs, SingleUnits, PlotUnits)

%% Process replay trials
if any(strcmp(TrialInfo.Perturbation,'OL-Template'))
    [OpenLoop] = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);
    ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotfigures', 1, 'whichunits', PlotUnits);
end

end
