%MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up
MySession = '/mnt/data/Processed/Behavior/PCX4/PCX4_20210721_r0_processed.mat'; % session path - leave empty to get browser pop up
LoadProcessedSession; % loads relevant variables

% any replays?
if any(strcmp(TrialInfo.Perturbation,'OL-Template'))
    [OpenLoop] = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);
else
    [OpenLoop] = [];
end

NumUnits = size(SingleUnits,2);
% sort units by tetrode - to match session viewer
clear foo
for i = 1:NumUnits
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id];
end
[~,SortedByTetrodes] = sort(foo(:,1));

PlotUnits = SortedByTetrodes; % plot all units
x = SortedByTetrodes(find(foo(SortedByTetrodes,1)==13));
y = SortedByTetrodes(find(foo(SortedByTetrodes,1)==14));
PlotUnits = [x; y];
PlotUnits = SortedByTetrodes(1:5);
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotfigures', 1, 'whichunits', PlotUnits);