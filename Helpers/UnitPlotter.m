function [x, AlignedFRs, BinOffset, RawSpikeCounts] = UnitPlotter(whichUnit, whichodor, whichTZ, AlignedSpikes, Events, TrialInfo, AlignTo, trialsdone, trialboxcolor)

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
whichUnit = abs(whichUnit);
perturbed = whichTZ<0; % hack to select perturbed or normal trials
whichTZ = abs(whichTZ);

thisUnitSpikes = AlignedSpikes(:,whichUnit);

whichTrials = intersect(find(TrialInfo.Odor==whichodor),find(TrialInfo.TargetZoneType==whichTZ));
if ~perturbed
    whichTrials = intersect(whichTrials, ...
        find(cellfun(@isempty, TrialInfo.Perturbation)));
else
    whichTrials = intersect(whichTrials, ...
        find(~cellfun(@isempty, TrialInfo.Perturbation)));
end

whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) TrialInfo.Duration(whichTrials)]; %#ok<AGROW>
whichTrials = sortrows(whichTrials,3);

AlignedFRs = []; RawSpikeCounts = []; BinOffset = NaN; x = 0;

if any(whichTrials)
    % Plot all events
    myEvents = Events(whichTrials(:,1),:);
    switch AlignTo
        case 1 % to trial ON
            Xlims = [-1.2 -1];
            Offset = 0*myEvents(:,1);
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
    
    if plotting
        % Plot Events
        EventPlotter(myEvents,trialsdone);
        % Plot TrialType
        TrialTypePlotter(whichTrials(:,2),whichodor,Xlims,trialsdone);
        
        % Plot Spikes
        for x = 1:size(whichTrials,1)
            % Plot Target Zone periods - adjust times if needed
            ZoneTimes = TrialInfo.InZone{whichTrials(x)} - Offset(x);
            InZonePlotter(ZoneTimes', x+trialsdone);
            
            % Plot Spikes
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
            
            % adjust spiketimes if needed
            thisTrialSpikes = thisTrialSpikes - Offset(x);
            PlotRaster(thisTrialSpikes,x+trialsdone,Plot_Colors('k'));
        end
    end
    
    % calculate PSTH

    BinOffset = Xlims(1)*1000;
    
    for TZ = 1:12
        thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1));
        Events2Align = Offset(find(whichTrials(:,2)==TZ),1);
        [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
        AlignedFRs(TZ,1:numel(myFR)) = myFR;
        RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
    end
end
%%
    function EventPlotter(myEvents,trialsdone)
        ticklength = 0.8;
        Y = trialsdone+((1:size(myEvents,1))' + repmat([-ticklength 0 0],size(myEvents,1),1) )';
        
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
        for j = 1:numel(U)
            y2 = y1 + numel(find(TrialList==U(j)));
            Y = [y1 y1 y2 y2];
            if perturbed
                fill(X,Y,trialboxcolor,'EdgeColor','k');
            else
                fill(X,Y,trialboxcolor,'EdgeColor','none');
            end
            y1 = y2;
        end
    end

end