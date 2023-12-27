function [x, FR, BinOffset] = PlotRandomTuningTrials(trialsdone, whichUnit, whichOdor, AlignedSpikes, Events, TrialSequence, AlignTo, LocationDuration, XLims, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;
psth = params.Results.psth;

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
if ~plotting % only for analysis
    plotevents = 0;
    plotspikes = 0;
end 
% colormap
%[whichcolor] = GetLocationColor(x);

whichUnit = abs(whichUnit);
thisUnitSpikes = AlignedSpikes(:,whichUnit);
% get trials of the particular odor
whichTrials = find(TrialSequence(:,2)==(whichOdor+1));
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

Xlims = [-1 2];
SpikesPSTH = [];

if plotevents
    line([-1 6], trialsdone + [0 0], 'Color', 'k');
    % Plot Events
    %EventPlotter(myEvents);
    % Plot TrialType
    TrialTypePlotter(whichTrials(:,1),whichOdor,[-1.2 -1],trialsdone);
end

for x = 1:size(whichTrials,1)
    if plotevents
        % Plot Odor On periods - adjust times if needed
        ZoneTimes = [0 LocationDuration];
        InZonePlotter(ZoneTimes', (x + trialsdone));
    end
    
    if plotspikes || ~plotting
        % Get Spikes
        thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
        
        % adjust spiketimes if needed
        thisTrialSpikes = thisTrialSpikes - Offset(x);
        
        % Plot Spikes
        if plotspikes
            PlotRaster(thisTrialSpikes,x + trialsdone,Plot_Colors('o'));
        end
        
        %SpikesPSTH = vertcat(SpikesPSTH, [thisTrialSpikes' (x)*ones(numel(thisTrialSpikes),1)]);
        SpikesPSTH = vertcat(SpikesPSTH, [thisTrialSpikes' (whichTrials(x,1))*ones(numel(thisTrialSpikes),1)]); %trying to see if same issue as before
    end
end

% thisodorSpikes = thisUnitSpikes(whichTrials);
% [AlignedFRs, RawSpikeCounts] = MakePSTH_v3(thisodorSpikes,Offset,BinOffset,'downsample',500);

% calculate PSTH
FR = [];
BinOffset = round(Xlims(1)*1000);
if psth
    % aggregate psth of all trials
    if ~isempty(SpikesPSTH)
        FR(1,:) = VerySimplePSTH(SpikesPSTH(:,1), x, BinOffset,'downsample',500);
    end
end

x = x + trialsdone;

    function InZonePlotter(TS, rowidx)
        foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.4);
        foo.EdgeColor = 'none';
        if ~isempty(TS)
            foo.Vertices = [ ...
                reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
                repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
            foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
        end
    end

    function TrialTypePlotter(TrialList,OdorType,Xlims,trialsdone)
        X = [Xlims(1) Xlims(2) Xlims(2) Xlims(1)];
        U = unique(TrialList);
        y1 = trialsdone;
        boxcolor(1,:) = Plot_Colors(['Odor',num2str(OdorType)]);
        boxcolor(2,:) = boxcolor-0.2;
        for j = 1:numel(U)
            y2 = y1 + numel(find(TrialList==U(j)));
            Y = [y1 y1 y2 y2];
            if trialsdone
                fill(X,Y,boxcolor(1,:),'EdgeColor','k','Linewidth',2);
            else
                fill(X,Y,boxcolor(1,:),'EdgeColor','none');
            end
            y1 = y2;
            boxcolor = flipud(boxcolor);
        end
    end

end