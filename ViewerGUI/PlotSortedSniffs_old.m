function [] = PlotSortedSniffs_old(whichUnit, whichOdor, SniffAlignedSpikes, TrialAlignedSniffs, TrialInfo, varargin)

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
whichUnit = abs(whichUnit);

Xlims = -[0.5 1.5];

% hack to prevent OL-Template trials to be considered as perturbed trials
f = find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'));
if ~isempty(f)
    for i = 1:numel(f)
        TrialInfo.Perturbation{f(i),1} = [];
    end
end

thisUnitSpikes = SniffAlignedSpikes(:,whichUnit);
whichodor = whichOdor;
% get the trials for that particular odor
allTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
    find(TrialInfo.Odor==whichodor));

% assemble a list of sniffs 
% [inh-start inh-end next-inh sniffID snifftype sniff-duration inh-duration Trial ID]

AllSniffs = [];
for s = 1:numel(allTrials) % every trial
    thisTrialSniffs = TrialAlignedSniffs{allTrials(s)}; 
    thisTrialSniffs(:,6) = thisTrialSniffs(:,3) - thisTrialSniffs(:,1);
    thisTrialSniffs(:,7) = thisTrialSniffs(:,2) - thisTrialSniffs(:,1);
    thisTrialSniffs(:,8) = s;
    AllSniffs = vertcat(AllSniffs, thisTrialSniffs);
end

% sort sniff List by Sniff Type, then sniff duration, then inh duration
AllSniffs = sortrows(AllSniffs,[5 6 7 8]);

if plotting
    
    if plotevents
        % Plot SniffType
        SniffTypePlotter(AllSniffs(:,5),whichodor,Xlims,0);
    end
    
    % Plot Spikes
    for x = 1:size(AllSniffs,1)
        if plotevents
            % Plot Target Zone periods - adjust times if needed
            ExhalationTimes = AllSniffs(x,[2 3]) - AllSniffs(x,1);
            SniffPlotter(ExhalationTimes', x);
        end
        
        if plotspikes
            
            if AllSniffs(x,4)<0
                whichtrial = AllSniffs(x,8) - 1;
                whichsniff = AllSniffs(x,4);
            else
                whichtrial = AllSniffs(x,8);
                whichsniff = AllSniffs(x,4);
            end
            
            if whichtrial
            % Plot Spikes
            thisTrialSpikes = thisUnitSpikes{whichtrial}{1};
            thisSniffSpikes = thisTrialSpikes(floor(thisTrialSpikes)==whichsniff) - whichsniff;
            
            % adjust spiketimes if needed
            PlotRaster(thisSniffSpikes,x,Plot_Colors('k'));
            end
        end
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