function [] = QuickCheckReplay(MySessionName)
SetupSmellocatorGlobals;
Paths = WhichComputer();
foo = regexp(MySessionName,'_','split');
AnimalName = foo{1};
MySession = fullfile(Paths.Grid.Behavior_processed,AnimalName,MySessionName);

% Load the relevant variables
load(MySession, 'Traces', 'TrialInfo', 'TTLs', 'ReplayTTLs', 'SingleUnits', 'PassiveReplayTraces', 'TargetZones');
PlotUnits = 8 + [1:8];
[OpenLoop] = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

% code to plot respiration traces
ReplayRespirationTracePlotter(OpenLoop.TemplateTraces, OpenLoop.ReplayTraces, PassiveReplayTraces);
set(gcf,'Position',[2615 27 700 969]);

ReplayTracePlotter(OpenLoop, TrialInfo, SingleUnits, TTLs, PassiveReplayTraces);

PlotUnits = [];
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
    'plotfigures', 1, 'plotephys', 1, 'whichunits', PlotUnits);
end