function [x,FR,BinOffset] = PlotSortedSniffs(whichUnit, whichOdor, TrialAlignedSpikes, TrialAlignedSniffs, TrialInfo, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('sortorder', 0, @(x) isnumeric(x)); % 0 - sniff duration, 1 - inhalation duration
params.addParameter('alignto', 1, @(x) isnumeric(x)); % 1 - inhalation start, 2 - inhalation end
params.addParameter('warptype', 0, @(x) isnumeric(x)); % 0 - no warp, 1 - by sniff duration, 2 - by inhalation duration
params.addParameter('BehaviorSampleRate', 500, @(x) isnumeric(x)); % 0 - sniff duration, 1 - inhalation duration

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;
psth = params.Results.psth;
sortby = params.Results.sortorder;
alignto = params.Results.alignto;
warptype = params.Results.warptype;
BehaviorSampRate = params.Results.BehaviorSampleRate;

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
whichUnit = abs(whichUnit);

Xlims = -[0.1 1.1];

% hack to prevent OL-Template trials to be considered as perturbed trials
f = find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'));
if ~isempty(f)
    for i = 1:numel(f)
        TrialInfo.Perturbation{f(i),1} = [];
    end
end

thisUnitSpikes = TrialAlignedSpikes(:,whichUnit);
whichodor = whichOdor;
% get the trials for that particular odor
allTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
    find(TrialInfo.Odor==whichodor));

% get halt trials if any
% also collect perturbation trials
perturbationTrials = intersect(find(strncmpi(TrialInfo.Perturbation(:,1),'Halt-Flip',9)), ...
    find(TrialInfo.Odor==whichodor));

% assemble a list of sniffs 
% [inh-start inh-end next-inh sniffID snifftype sniff-duration inh-duration Trial ID]

AllSniffs = [];
for s = 1:numel(allTrials) % every trial
    thisTrialSniffs = TrialAlignedSniffs{allTrials(s)}; 
    thisTrialSniffs(:,15) = thisTrialSniffs(:,3) - thisTrialSniffs(:,1); % col 15 = sniff duration
    thisTrialSniffs(:,16) = thisTrialSniffs(:,2) - thisTrialSniffs(:,1); % col 16 = inhalation duration
    thisTrialSniffs(:,17) = allTrials(s); % col 17 = trial ID
    AllSniffs = vertcat(AllSniffs, thisTrialSniffs);
end

% add perturbation sniffs if any
if ~isempty(perturbationTrials)
    for s = 1:numel(perturbationTrials)
        thisTrialSniffs = TrialAlignedSniffs{perturbationTrials(s)}; 
        haltperiod = TrialInfo.Perturbation{perturbationTrials(s),2}(1:2)./BehaviorSampRate;
        haltsniffs = intersect(find(thisTrialSniffs(:,1)>=haltperiod(1)),find(thisTrialSniffs(:,1)<=haltperiod(2)));
        thisTrialSniffs = thisTrialSniffs(haltsniffs,:);
        thisTrialSniffs(:,15) = thisTrialSniffs(:,3) - thisTrialSniffs(:,1);
        thisTrialSniffs(:,16) = thisTrialSniffs(:,2) - thisTrialSniffs(:,1);
        thisTrialSniffs(:,17) = perturbationTrials(s);
        % change snifftype 
        thisTrialSniffs(:,5) = 3;
        AllSniffs = vertcat(AllSniffs, thisTrialSniffs);
    end
end

switch sortby
    case 0
        % sort sniff List by Sniff Type, then sniff duration, then inh duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[5 15 16 17]);
    case 1
        % sort sniff List by Sniff Type, then inh duration, then sniff duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[5 16 15 17]);
    case 2
        AllSniffs(find(AllSniffs(:,5)==0),5) = 1; 
        AllSniffs(find(AllSniffs(:,5)==1),5) = 1;
        AllSniffs(find(AllSniffs(:,5)==2),5) = 1; 
        % sort by sniff type, then odor location (col 6), then sniff duration, then inh duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[5 6 16 15 17]);
end

SpikesPlot = []; SpikesPSTH = [];
if plotting
    
    % Plot Spikes
    for x = 1:size(AllSniffs,1)
        if plotevents
            % Plot Target Zone periods - adjust times if needed
            ExhalationTimes = AllSniffs(x,[8 9 2 3 12 13]) - AllSniffs(x,alignto);
            ExhalationTimes = reshape(ExhalationTimes,2,3)';
            if warptype
                ExhalationTimes = ExhalationTimes * (mean(AllSniffs(:,14+warptype))/AllSniffs(x,14+warptype));
            end
            SniffPlotter(ExhalationTimes', x);
        end
        
        if plotspikes
              
              whichtrial = AllSniffs(x,end);
              %whichsniff = AllSniffs(x,4);
              
              thisTrialSpikes = thisUnitSpikes{whichtrial}{1} - AllSniffs(x,alignto);
              if warptype
                  thisTrialSpikes = thisTrialSpikes * (mean(AllSniffs(:,14+warptype))/AllSniffs(x,14+warptype));
              end
              
              SpikesPlot = vertcat(SpikesPlot, [thisTrialSpikes' x*ones(numel(thisTrialSpikes),1)]);

              % for plotting PSTH
              switch alignto
                  case 1 % inhalation start
                      thisTrialSpikes(:,thisTrialSpikes>AllSniffs(x,15)) = [];
                  case 2 % inhalation end
                      thisTrialSpikes(:,thisTrialSpikes>(AllSniffs(x,15)-AllSniffs(x,16))) = [];
              end
              SpikesPSTH = vertcat(SpikesPSTH, [thisTrialSpikes' x*ones(numel(thisTrialSpikes),1)]);

        end
    end
end

BinOffset = -1000;
if plotspikes
    plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize',0.5);
end

x = size(AllSniffs,1);
if psth
    for i = -1:1:3
        whichsniffs = find(AllSniffs(:,5)==i);
        if ~isempty(whichsniffs)
            %FR{i+2} = MakePSTH_v4(SpikesPlot(find(ismember(SpikesPlot(:,2),whichsniffs)),1),numel(whichsniffs),BinOffset,'downsample',500,'kernelsize',20);
            FR{i+2} = MakeSniffTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),whichsniffs)),1),...
                AllSniffs(whichsniffs,15),...
                BinOffset,'downsample',500,'kernelsize',20);
        end
    end
else
    FR = [];
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

    function SniffPlotter(TS, rowidx)
        foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.4);
        foo.EdgeColor = 'none';
        if ~isempty(TS)
            foo.Vertices = [ ...
                reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
                repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
            foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
        end
    end

    function SniffTypePlotter(TrialList,OdorType,Xlims,trialsdone)
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