function [] = QuickCheckReplay(MySessionName)
SetupSmellocatorGlobals;
Paths = WhichComputer();
foo = regexp(MySessionName,'_','split');
AnimalName = foo{1};
MySession = fullfile(Paths.Grid.Behavior_processed,AnimalName,MySessionName);
global pdfPosition
pdfPosition = [0.6336    0.0574    0.3086    0.8028];
global MyPDFname

MyPDFname = '/mnt/data/Sorted/Q9/2022-11-19_17-14-37/ClusterMaps/ReplayResponses.pdf';

% Load the relevant variables
load(MySession, 'Traces', 'TrialInfo', 'TTLs', 'ReplayTTLs', 'SingleUnits', 'PassiveReplayTraces', 'TargetZones');
PlotUnits = 8 + [1:8];
PlotUnits = [];
if any(strcmp(TrialInfo.Perturbation(:,1),'OL-Replay'))
    [OpenLoop] = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

    % code to plot respiration traces
    %ReplayRespirationTracePlotter(OpenLoop.TemplateTraces, OpenLoop.ReplayTraces, PassiveReplayTraces);
    ReplaySelectTracePlotter(OpenLoop.TemplateTraces, OpenLoop.ReplayTraces, PassiveReplayTraces,'Sniffs');
    set(gcf,'Position',[2615 27 700 969]);

    %ReplayTracePlotter(OpenLoop, TrialInfo, SingleUnits, TTLs, PassiveReplayTraces);

    ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
    'plotfigures', 1, 'savepdfs', 1, 'plotephys', 1, 'whichunits', PlotUnits);
else
    [OpenLoop] = ExtractReplayTrials_v2(Traces, TrialInfo, TTLs, ReplayTTLs, PassiveReplayTraces);

    ProcessOpenLoopTrials_v2(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
    'plotfigures', 1, 'plotephys', 1, 'whichunits', PlotUnits);
end

end