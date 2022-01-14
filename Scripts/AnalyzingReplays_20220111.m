MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up
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
PlotUnits = SortedByTetrodes(find(foo(SortedByTetrodes,1)==7));
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotfigures', 1, 'whichunits', PlotUnits);