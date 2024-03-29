function [x, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = PlotFullSession(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, ZoneTimesIn, AlignTo, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('sniffaligned', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('sniffscalar', 3, @(x) isnumeric(x));
params.addParameter('trialfilter', 1, @(x) isnumeric(x)); % 1 - plot all trials, 2 - plot only halt matched TZs, 3 - only use the hlat templates
params.addParameter('poolTZs', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;
psth = params.Results.psth;
sniffaligned = params.Results.sniffaligned;
sniffscalar = params.Results.sniffscalar;
trialfilter = params.Results.trialfilter;
poolTZs = params.Results.poolTZs; % group trials of different TZs together when calculating PSTH

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
whichUnit = abs(whichUnit);
thisUnitSpikes = AlignedSpikes(:,whichUnit);
whichodor = whichOdor;

if trialfilter == 3 % only use template trials
    whichTrials = intersect(find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template')), ...
        find(TrialInfo.Odor==whichodor));
    
    % hack to prevent OL-Template trials to be considered as perturbed trials
    f = find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'));
    if ~isempty(f)
        for i = 1:numel(f)
            TrialInfo.Perturbation{f(i),1} = [];
        end
    end
    
else
    % hack to prevent OL-Template trials to be considered as perturbed trials
    f = find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'));
    if ~isempty(f)
        for i = 1:numel(f)
            TrialInfo.Perturbation{f(i),1} = [];
        end
    end
    
    % get the trial sorting order
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
        find(TrialInfo.Odor==whichodor));
    
end
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) TrialInfo.Duration(whichTrials)];
if isfield(TrialInfo,'SessionSplits')
    whichTrials(:,4) = whichTrials(:,1)>TrialInfo.SessionSplits(1,2);
else
    whichTrials(:,4) = 0;
end
whichTrials = sortrows(whichTrials,[4,2,3]); % to keep individual closed loop sessions separate 

% also collect perturbation trials
perturbationTrials = intersect(find(~cellfun(@isempty, TrialInfo.Perturbation)), ...
    find(TrialInfo.Odor==whichodor));
perturbationTrials = intersect(find(~strcmp(TrialInfo.Perturbation(:,1),'OL-Replay')), ...
    perturbationTrials);
perturbationTrials = [perturbationTrials TrialInfo.TargetZoneType(perturbationTrials) TrialInfo.Duration(perturbationTrials)]; %#ok<AGROW>
perturbationTrials(:,4) = 0; % hack for dealing with separation of concatenated closed-loop sessions
perturbationTrials = sortrows(perturbationTrials,[2,3]);

% for offsets - sort by offset type
if any(perturbationTrials)
    if any(strcmp(TrialInfo.Perturbation(:,1),'Offset-II-Template')) || ...
            any(strcmp(TrialInfo.Perturbation(:,1),'Offset-II'))
        offsetParams = cell2mat(TrialInfo.Perturbation(perturbationTrials(:,1),2));
        [~,sortidx] = sort(offsetParams(:,3));
        perturbationTrials = perturbationTrials(sortidx,:);
        offsetParams = offsetParams(sortidx,:);
        offsettypes = unique(offsetParams(:,3));
        for k = 1:numel(offsettypes)
            f = find(offsetParams(:,3)==offsettypes(k));
            perturbationTrials(f,2) = perturbationTrials(f,2) + 0.1*k;
        end
    end
    perturbationTZs = unique(perturbationTrials(:,2));
else
    perturbationTZs = [];
end

if trialfilter == 2 % only keep close loop trials that match the halted TZs
    f = find(~ismember(whichTrials(:,2),perturbationTZs));
    whichTrials(f,:) = [];
end

allTrials = vertcat(whichTrials, perturbationTrials);

% Plot all events
myEvents = Events(allTrials(:,1),:);
switch AlignTo
    case 1 % to trial ON
        Xlims = [-1.2 -1];
        Offset = 0*myEvents(:,1);
    case 2 % odor ON
        if ~sniffaligned
            odorON = myEvents(:,1);
        else
            odorON = floor(myEvents(:,1));
        end
        myEvents(:,1) = 0; % replace odorON with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - odorON;
        Xlims = [-1.2 -1];
        Offset = odorON;
    case 3 % trial OFF
        if ~sniffaligned
            TrialOFF = myEvents(:,3);
        else
            TrialOFF = floor(myEvents(:,3));
        end
        myEvents(:,3) = 0; % replace TrialOFF with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - TrialOFF;
        Xlims = [-1.2 -1] - 4;
        Offset = TrialOFF;
    case 4 % reward
        if ~sniffaligned
            Reward = myEvents(:,3);
        else
            Reward = myEvents(:,3);
        end
        myEvents(:,2) = 0; % replace Reward with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - Reward;
        Xlims = [-1.2 -1] - 4;
        Offset = Reward;
    case 5 % first TZ entry with stay > 100ms
        if ~sniffaligned
            Offset = myEvents(:,4);
        else
            Offset = floor(myEvents(:,4));
        end
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1] - 1;
    case 6 % perturbation start
        if ~sniffaligned
            Offset = myEvents(:,5);
        else
            Offset = floor(myEvents(:,5));
        end
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1];
end
if sniffaligned
    Xlims = sniffscalar*Xlims;
end

lineplotted = 0;
if plotting
    
    if plotevents
        % Plot Events
        EventPlotter(myEvents);
        % Plot TrialType
        TrialTypePlotter(whichTrials(:,2),whichodor,Xlims,0);
    end
    
    % Plot Spikes
    for x = 1:size(whichTrials,1)
        if plotevents
            % Plot Target Zone periods - adjust times if needed
            ZoneTimes = ZoneTimesIn{whichTrials(x,1)} - Offset(x);
            InZonePlotter(ZoneTimes', x);
            
            if ~lineplotted && whichTrials(x,4) == 1
                line([-1.5 6], [x x], 'Color', 'k');
                lineplotted = lineplotted + 1;
            end
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
BinOffset = round(Xlims(1)*1000);
if psth
    if ~poolTZs
        for TZ = 1:12
            thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1));
            Events2Align = Offset(find(whichTrials(:,2)==TZ),1);
            [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
            AlignedFRs(TZ,1:numel(myFR)) = myFR;
            RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
        end
    else
        TZ = 1;
        thisTZspikes = thisUnitSpikes(whichTrials(:,1));
        Events2Align = Offset(find(whichTrials(:,2)),1);
        [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
        AlignedFRs(TZ,1:numel(myFR)) = myFR;
        RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH; 
    end
    entries_done = TZ;
end

if ~isempty(perturbationTrials)
    if plotting
        if plotevents
            TrialTypePlotter(perturbationTrials(:,2),whichodor,Xlims,x);
        end
        
        
        for y = 1:size(perturbationTrials,1)
            if plotevents
                % Plot Target Zone periods - adjust times if needed
                ZoneTimes = ZoneTimesIn{perturbationTrials(y)} - Offset(x+y);
                InZonePlotter(ZoneTimes', y+x);
            end
            if plotspikes
                % Plot Spikes
                thisTrialSpikes = thisUnitSpikes{perturbationTrials(y,1)}{1};
                % adjust spiketimes if needed
                thisTrialSpikes = thisTrialSpikes - Offset(x+y);
                PlotRaster(thisTrialSpikes,x+y,Plot_Colors('Paletton',[1 2]));
            end
        end
    end
end

% calculate PSTH
AlignedPerturbationFRs = [];
if psth
    if ~poolTZs
        for TZ = 1:12
            thisTZspikes = thisUnitSpikes(perturbationTrials(find(perturbationTrials(:,2)==TZ),1));
            Events2Align = Offset(x+find(perturbationTrials(:,2)==TZ),1);
            [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
            AlignedPerturbationFRs(TZ,1:numel(myFR)) = myFR;
            RawSpikeCounts(TZ+entries_done,1:numel(myPSTH)) = myPSTH;
        end
    else
        TZ = 1;
        thisTZspikes = thisUnitSpikes(perturbationTrials(:,1));
        Events2Align = Offset(x+find(perturbationTrials(:,2)),1);
        [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
        AlignedPerturbationFRs(TZ,1:numel(myFR)) = myFR;
        RawSpikeCounts(TZ+entries_done,1:numel(myPSTH)) = myPSTH;
    end
end

x = size(allTrials,1);

%%
    function EventPlotter(myEvents)
        ticklength = 0.8;
        Y = ((1:size(myEvents,1))' + repmat([-ticklength 0 0],size(myEvents,1),1) )';
        
        % Odor ON
        X = [repmat(myEvents(:,1),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('r'),'Linewidth',2);
        
        % Trial OFF
        X = [repmat(myEvents(:,3),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('r'),'Linewidth',2);
        
        % Rewards
        X = [repmat(myEvents(:,2),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('b'),'Linewidth',2);
    end

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
                fill(X,Y,boxcolor(1,:),'EdgeColor','k');
            else
                fill(X,Y,boxcolor(1,:),'EdgeColor','none');
            end
            y1 = y2;
            boxcolor = flipud(boxcolor);
        end
    end
end