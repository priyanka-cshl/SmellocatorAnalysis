function [x, ActiveFRs, PassiveFRs, RawSpikeCounts] = AddReplay2FullSession(trialsdone, whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo, SortTrials, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('poolTZs', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;
psth = params.Results.psth;
poolTZs = params.Results.poolTZs; % group trials of different TZs together when calculating PSTH

thisUnitSpikes = AlignedSpikes(:,whichUnit);
whichodor = whichOdor;
% get the trial sorting order
whichTrials = find(TrialInfo.Odor==whichodor); % both active and passive replays
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) ...
    TrialInfo.Duration(whichTrials) (TrialInfo.TrialID(whichTrials)<0)']; %#ok<AGROW>

if isfield(TrialInfo,'SessionId')
    % change column 4 to be 0 - first active replay, 
    % 1 - passive replay 1, 2 -
    % passive replay 2 etc
    whichTrials(:,4) = TrialInfo.SessionId(whichTrials(:,1));
end

% Sort trials - first by active and passive replay
if SortTrials
    whichTrials = sortrows(whichTrials,[4,2,3]); % to keep individual replay sessions separate
end

% Plot all events
myEvents = Events(whichTrials(:,1),:);
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
end

if plotevents
    EventPlotter(myEvents,trialsdone);
    
    % Plot TrialType
    TrialTypePlotter(whichTrials(:,2),whichodor,Xlims,trialsdone);
end

% Plot Spikes
lineplotted = 0;
for x = 1:size(whichTrials,1)
    if plotevents
        % Plot Target Zone periods - adjust times if needed
        ZoneTimes = TrialInfo.InZone{whichTrials(x)} - Offset(x);
        InZonePlotter(ZoneTimes', x+trialsdone);
        if ~lineplotted
            line([-1.5 6], trialsdone - 1 + [x x], 'Color', 'k');
            lineplotted = lineplotted + 1;
        elseif lineplotted == 1 && whichTrials(x,4) == 1 % likely first passive replay round
            line([-1.5 6], trialsdone -1 + [x x], 'Color', 'k');
            lineplotted = lineplotted + 1;
        elseif lineplotted == 2 && whichTrials(x,4) == 2 % likely second passive replay round
            line([-1.5 6], trialsdone -1 + [x x], 'Color', 'k');
            lineplotted = lineplotted + 1;
        end
        
    end
    
    if plotspikes
        % Plot Spikes
        thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
        % adjust spiketimes if needed
        thisTrialSpikes = thisTrialSpikes - Offset(x);
        
        if TrialInfo.TrialID(whichTrials(x))>0
            PlotRaster(thisTrialSpikes,x+trialsdone,Plot_Colors('r'));
        else
            PlotRaster(thisTrialSpikes,x+trialsdone,Plot_Colors('t'));
        end
    end
end

%%
x = size(whichTrials,1);
% calculate PSTH
ActiveFRs = []; PassiveFRs = []; RawSpikeCounts = []; 
BinOffset = round(Xlims(1)*1000);

if psth
    % active replays
    if ~poolTZs
        for TZ = 1:12
            whichreplaytrials = intersect(find(whichTrials(:,2)==TZ),find(whichTrials(:,4)==0));
            thisTZspikes = thisUnitSpikes(whichTrials(whichreplaytrials,1));
            Events2Align = Offset(whichreplaytrials,1);
            [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
            ActiveFRs(TZ,1:numel(myFR)) = myFR;
            RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
        end
    else
        TZ = 1;
        whichreplaytrials = find(whichTrials(:,4)==0); % all TZs
        thisTZspikes = thisUnitSpikes(whichTrials(whichreplaytrials,1));
        Events2Align = Offset(whichreplaytrials,1);
        [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
        ActiveFRs(TZ,1:numel(myFR)) = myFR;
        RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH; 
    end
    entries_done = TZ;
    
    % passive replays
    if ~poolTZs
        for TZ = 1:12
            whichreplaytrials = intersect(find(whichTrials(:,2)==TZ),find(whichTrials(:,4)==1));
            thisTZspikes = thisUnitSpikes(whichTrials(whichreplaytrials,1));
            Events2Align = Offset(whichreplaytrials,1);
            [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
            PassiveFRs(TZ,1:numel(myFR)) = myFR;
            RawSpikeCounts(TZ+entries_done,1:numel(myPSTH)) = myPSTH;
        end
    else
        TZ = 1;
        whichreplaytrials = find(whichTrials(:,4)==1); % all TZs
        thisTZspikes = thisUnitSpikes(whichTrials(whichreplaytrials,1));
        Events2Align = Offset(whichreplaytrials,1);
        [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
        PassiveFRs(TZ,1:numel(myFR)) = myFR;
        RawSpikeCounts(TZ+entries_done,1:numel(myPSTH)) = myPSTH; 
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
        boxcolor(1,:) = Plot_Colors(['Odor',num2str(OdorType)]);
        boxcolor(2,:) = boxcolor-0.2;
        for j = 1:numel(U)
            y2 = y1 + numel(find(TrialList==U(j)));
            Y = [y1 y1 y2 y2] + trialsdone;
            fill(X,Y,boxcolor(1,:),'EdgeColor','none');
            y1 = y2;
            boxcolor = flipud(boxcolor);
        end
    end
end