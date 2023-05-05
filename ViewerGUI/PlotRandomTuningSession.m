function [x, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = PlotRandomTuningSession(whichUnit, whichOdor, AlignedSpikes, Events, TrialSequence, AlignTo, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;

% colormap
cmap = colormap(
plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
whichUnit = abs(whichUnit);

thisUnitSpikes = AlignedSpikes(:,whichUnit);
whichodor = whichOdor;
% get trials of the particular odor
whichTrials = find(TrialSequence(:,2)==whichodor);

% create a Events matrix - [LocationStart LocationEnd Location]
Offset(1:whichTrials) = 0; % for alignment

if plotting
    
    % Plot Spikes
    for x = 1:size(whichTrials,1)
        if plotevents
            % Plot Target Zone periods - adjust times if needed
            LocationTimes = [Events - Offset(x) TrialSequence(whichTrials(x),3:end)'];
            LocationPlotter(LocationTimes, x);
        end
        if plotspikes
            % Plot Spikes
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
            
            % adjust spiketimes if needed
            thisTrialSpikes = thisTrialSpikes - Offset(x);
            PlotRaster(thisTrialSpikes,x,Plot_Colors('k'));
        end
    end
end

x = size(whichTrials,1);
% calculate PSTH
AlignedFRs = []; RawSpikeCounts = [];
BinOffset = Xlims(1)*1000;
x = size(allTrials,1);

%%


    function LocationPlotter(TS, rowidx)
        foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.4);
        foo.EdgeColor = 'none';
        if ~isempty(TS)
            foo.Vertices = [ ...
                reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
                repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
            foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
        end
    end

end