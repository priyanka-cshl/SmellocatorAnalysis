function [AlignedFRs, BinOffset, RawSpikeCounts] = PlotRandomTuningSession(whichUnit, whichOdor, AlignedSpikes, Events, TrialSequence, AlignTo, XLims, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis

% colormap
%[whichcolor] = GetLocationColor(x);

whichUnit = abs(whichUnit);
thisUnitSpikes = AlignedSpikes(:,whichUnit);
% get trials of the particular odor
whichTrials = find(TrialSequence(:,2)==whichOdor);
BinOffset = XLims(1)*1000;

LMat = TrialSequence(whichTrials,3:end);
if AlignTo == 1 % Trial ON
    Offset(1:numel(whichTrials)) = 0;
end
if AlignTo == 2 % odor ON
    Offset(1:numel(whichTrials)) = 0;
end
if AlignTo>2 % align to a specific location
    whichLocation = 1000-AlignTo;
    for x = 1:size(whichTrials,1)
        idx = find(LMat(x,:)==whichLocation);
        Offset(x) = Events.LocationShifts(idx,1);
        if Offset(x)<0
            Offset(x) = Events.Odor(1);
        end
        LMat(x,:) = circshift(LMat(x,:),-idx + ceil(size(LMat,2)/2));
    end
end

if plotting
    if plotevents
        image(LMat+250);
    end
    
    if plotspikes
        for x = 1:size(whichTrials,1)
            % Plot Spikes
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
            
            % adjust spiketimes if needed
            thisTrialSpikes = thisTrialSpikes - Offset(x);
            PlotRaster(thisTrialSpikes,x,Plot_Colors('k'));
        end
    end
end

thisodorSpikes = thisUnitSpikes(whichTrials);
[AlignedFRs, RawSpikeCounts] = MakePSTH_v3(thisodorSpikes,Offset,BinOffset,'downsample',500);

end