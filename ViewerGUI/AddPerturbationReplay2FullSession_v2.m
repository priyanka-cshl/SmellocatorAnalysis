function [x, PassiveReplayFRs, PerturbationReplayFRs, BinOffset] = ...
    AddPerturbationReplay2FullSession_v2(trialsdone, whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, ZoneTimesIn, AlignTo, SortTrials, varargin)

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

% get the trial sorting order
whichTrials = intersect(find(TrialInfo.Odor==whichodor), find(~TrialInfo.Perturbed));
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) ...
    TrialInfo.Duration(whichTrials)]; %#ok<AGROW>

% Sort trials by target zone type
if SortTrials
    for tz = 1:12
        whichTrials(whichTrials(:,2)==tz,:) = sortrows(whichTrials(whichTrials(:,2)==tz,:),3);
    end
end

% also collect perturbation trials
perturbationTrials = intersect(find(TrialInfo.Odor==whichodor), find(TrialInfo.Perturbed));
perturbationTrials = [perturbationTrials TrialInfo.TargetZoneType(perturbationTrials) TrialInfo.Duration(perturbationTrials)]; %#ok<AGROW>
perturbationTrials = sortrows(perturbationTrials,2);
for tz = 1:12
    perturbationTrials(perturbationTrials(:,2)==tz,:) = sortrows(perturbationTrials(perturbationTrials(:,2)==tz,:),3);
end

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

allTrials = vertcat(perturbationTrials, whichTrials);

% Plot all events
myEvents = Events(allTrials(:,1),:);
switch AlignTo
    case 1 % to trial ON
        Xlims = [-1.2 -1];
        Offset = 0*myEvents(:,4);
    case 2 % odor ON
        odorON = myEvents(:,1);
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
        TrialOFF = myEvents(:,3);
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
        Reward = myEvents(:,3);
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
        Offset = myEvents(:,4);
        if ~sniffaligned
            Offset = myEvents(:,4);
        else
            Offset = floor(myEvents(:,4));
        end
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1] - 1;
    case 6 % perturbation start
        Offset = myEvents(:,5);
        if ~sniffaligned
            Offset = myEvents(:,5);
        else
            Offset = floor(myEvents(:,5));
        end
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1];
end

if plotting
    if plotevents
        line([-1 6], trialsdone + [0 0], 'Color', 'k');
        
        EventPlotter(myEvents,trialsdone);
        
        % Plot TrialType
        TrialTypePlotter(perturbationTrials(:,2),-whichodor,Xlims,trialsdone);
        
        TrialTypePlotter(whichTrials(:,2),whichodor,Xlims,trialsdone+size(perturbationTrials,1));
    end
    
    if sniffaligned
        Xlims = sniffscalar*Xlims;
    end
    
    SpikesPSTH = [];
    
    % Plot Spikes
    for x = 1:size(allTrials,1)
        
        if plotevents
            % Plot Target Zone periods - adjust times if needed
            %ZoneTimes = TrialInfo.InZone{allTrials(x)} - Offset(x);
            ZoneTimes = ZoneTimesIn{allTrials(x)} - Offset(x);
            InZonePlotter(ZoneTimes', x+trialsdone);
        end
        
        if plotspikes
            % Plot Spikes
            thisTrialSpikes = thisUnitSpikes{allTrials(x,1)}{1};
            % adjust spiketimes if needed
            thisTrialSpikes = thisTrialSpikes - Offset(x);
            
            if x>size(perturbationTrials,1) % passive replays
                PlotRaster(thisTrialSpikes,x+trialsdone,Plot_Colors('b'));
            else % replayed perturbations
                PlotRaster(thisTrialSpikes,x+trialsdone,Plot_Colors('r'));
            end
            
            %SpikesPSTH = vertcat(SpikesPSTH, [thisTrialSpikes x*ones(numel(thisTrialSpikes),1)]);
            SpikesPSTH = vertcat(SpikesPSTH, [thisTrialSpikes (allTrials(x,1))*ones(numel(thisTrialSpikes),1)]);
        end
    end
end

% with current PG code I can't extract for analysis 
if plotting == 0
    SpikesPSTH = [];
for x = 1:size(allTrials,1)
    thisTrialSpikes = thisUnitSpikes{allTrials(x,1)}{1};
    % adjust spiketimes if needed
    thisTrialSpikes = thisTrialSpikes - Offset(x);
    SpikesPSTH = vertcat(SpikesPSTH, [thisTrialSpikes (allTrials(x,1))*ones(numel(thisTrialSpikes),1)]);
end 
end

PassiveReplayFRs = [];
PerturbationReplayFRs = [];
BinOffset = round(Xlims(1)*1000);

if psth && ~isempty(SpikesPSTH)
    % passive replays
    if poolTZs
        % passive replays
        PassiveReplayFRs{1} = MakeTrialTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),whichTrials(:,1))),1),...
            whichTrials(:,3),...
            BinOffset,'downsample',500);
        % perturbation replays
        PerturbationReplayFRs{1} = MakeTrialTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),perturbationTrials(:,1))),1),...
            perturbationTrials(:,3),...
            BinOffset,'downsample',500);
    else
        % TZ specific FRs
        for tz = 1:12
            % passive replays
            thisTZreplays = find(whichTrials(:,2)==tz);
            if ~isempty(find(ismember(SpikesPSTH(:,2),whichTrials(thisTZreplays,1))))
                PassiveReplayFRs{tz} = MakeTrialTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),whichTrials(thisTZreplays,1))),1),...
                    whichTrials(thisTZreplays,3),...
                    BinOffset,'downsample',500);
            end
            
            % perturbation replays
            thisPTreplays = find(perturbationTrials(:,2)==tz);
            if ~isempty(find(ismember(SpikesPSTH(:,2),perturbationTrials(thisPTreplays,1))))
                PerturbationReplayFRs{tz} = MakeTrialTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),perturbationTrials(thisPTreplays,1))),1),...
                    perturbationTrials(thisPTreplays,3),...
                    BinOffset,'downsample',500);
            end
            
        end
    end
end

%%
    function EventPlotter(myEvents,trialsdone)
        ticklength = 0.8;
        Y = (trialsdone+(1:size(myEvents,1))' + repmat([-ticklength 0 0],size(myEvents,1),1) )';
        
        % Odor ON
        X = [repmat(myEvents(:,1),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('pl'),'Linewidth',2);
        
        % Trial OFF
        X = [repmat(myEvents(:,3),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('pl'),'Linewidth',2);
        
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
        y1 = 0;
        boxcolor(1,:) = Plot_Colors(['Odor',num2str(abs(OdorType))]);
        boxcolor(2,:) = boxcolor-0.2;
        for j = 1:numel(U)
            y2 = y1 + numel(find(TrialList==U(j)));
            Y = [y1 y1 y2 y2] + trialsdone;
            if OdorType<0
                fill(X,Y,boxcolor(1,:),'EdgeColor','none');
            else
                fill(X,Y,boxcolor(1,:),'EdgeColor','k');
            end
            y1 = y2;
            boxcolor = flipud(boxcolor);
        end
    end
end