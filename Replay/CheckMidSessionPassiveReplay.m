function [] = CheckMidSessionPassiveReplay(MySession)

% Load the relevant variables
load(MySession, 'Traces', 'TrialInfo', 'TTLs', 'ReplayTTLs', 'SingleUnits');
[OpenLoop] = ExtractReplayTrials(Traces{1},TrialInfo{1},TTLs,ReplayTTLs);
PlotUnits = [1:8];
ProcessOpenLoopTrials(OpenLoop, TrialInfo{1}, SingleUnits, TTLs, ...
    'plotfigures', 1, 'plotephys', 1, 'whichunits', PlotUnits);

end