function [x] = PlotFullSession(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo)


thisUnitSpikes = AlignedSpikes(:,whichUnit);
whichodor = whichOdor;
% get the trial sorting order
whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
    find(TrialInfo.Odor==whichodor));
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) TrialInfo.Duration(whichTrials)]; %#ok<AGROW>
whichTrials = sortrows(whichTrials,2);

for tz = 1:12
    whichTrials(whichTrials(:,2)==tz,:) = sortrows(whichTrials(whichTrials(:,2)==tz,:),3);
end

% also collect perturbation trials
perturbationTrials = intersect(find(~cellfun(@isempty, TrialInfo.Perturbation)), ...
    find(TrialInfo.Odor==whichodor));
perturbationTrials = [perturbationTrials TrialInfo.TargetZoneType(perturbationTrials) TrialInfo.Duration(perturbationTrials)]; %#ok<AGROW>
perturbationTrials = sortrows(perturbationTrials,2);
for tz = 1:12
    perturbationTrials(perturbationTrials(:,2)==tz,:) = sortrows(perturbationTrials(perturbationTrials(:,2)==tz,:),3);
end

allTrials = vertcat(whichTrials, perturbationTrials);

% Plot all events
myEvents = Events(allTrials(:,1),:);
switch AlignTo
    case 1 % to trial ON
        Xlims = [-1.2 -1];
        Offset = 0*myEvents(:,4);
    case 2 % odor ON
        odorON = myEvents(:,1);
        myEvents(:,1) = 0; % replace odorON with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - odorON;
        Xlims = [-1.2 -1];
        Offset = odorON;
    case 3 % trial OFF
        TrialOFF = myEvents(:,3);
        myEvents(:,3) = 0; % replace TrialOFF with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - TrialOFF;
        Xlims = [-1.2 -1] - 4;
        Offset = TrialOFF;
    case 4 % reward
        Reward = myEvents(:,3);
        myEvents(:,2) = 0; % replace Reward with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - Reward;
        Xlims = [-1.2 -1] - 4;
        Offset = Reward;
    case 5 % first TZ entry with stay > 100ms
        Offset = myEvents(:,4);
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1] - 1;
    case 6 % perturbation start
        Offset = myEvents(:,5);
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1];
end

EventPlotter(myEvents);

% Plot TrialType
TrialTypePlotter(whichTrials(:,2),whichodor,Xlims,0);

% Plot Spikes
for x = 1:size(whichTrials,1)
    
    % Plot Target Zone periods - adjust times if needed
    ZoneTimes = TrialInfo.InZone{whichTrials(x)} - Offset(x);
    InZonePlotter(ZoneTimes', x);
    
    % Plot Spikes
    thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
    % adjust spiketimes if needed
    thisTrialSpikes = thisTrialSpikes - Offset(x);
    PlotRaster(thisTrialSpikes,x,Plot_Colors('k'));
end

if ~isempty(perturbationTrials)
TrialTypePlotter(perturbationTrials(:,2),whichodor,Xlims,x);
    for y = 1:size(perturbationTrials,1)
        
        % Plot Target Zone periods - adjust times if needed
        ZoneTimes = TrialInfo.InZone{perturbationTrials(y)} - Offset(x+y);
        InZonePlotter(ZoneTimes', y+x);
        
        % Plot Spikes
        thisTrialSpikes = thisUnitSpikes{perturbationTrials(y,1)}{1};
        % adjust spiketimes if needed
        thisTrialSpikes = thisTrialSpikes - Offset(x+y);
        PlotRaster(thisTrialSpikes,x+y,Plot_Colors('Paletton',[1 2]));
    end
end

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