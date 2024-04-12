function [] = QuickCheckReplay(MySessionName)
Paths = WhichComputer();
foo = regexp(MySessionName,'_','split');
AnimalName = foo{1};
MySession = fullfile(Paths.Grid.Behavior_processed,AnimalName,MySessionName);

% Load the relevant variables
load(MySession, 'Traces', 'TrialInfo', 'TTLs', 'ReplayTTLs', 'SingleUnits', 'PassiveReplayTraces');
PlotUnits = 8 + [1:8];
[OpenLoop] = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

ReplayTracePlotter(OpenLoop, TrialInfo, SingleUnits, TTLs, PassiveReplayTraces);

ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
    'plotfigures', 1, 'plotephys', 1, 'whichunits', PlotUnits);

end