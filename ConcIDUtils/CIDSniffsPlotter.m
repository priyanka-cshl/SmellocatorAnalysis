
[StimSettings, TTLs, SingleUnits, AllSpikes, TrialWiseSniffs, SniffsPlot] = CIDResponsePrepper(myKsDir);

% for saving figures to pdf
savefigs = 0;
if savefigs
    FigPath = strrep(myKsDir,'/mnt/data/Sorted','/home/priyanka/Desktop/cid');
    if ~exist("FigPath",'dir')
        mkdir(FigPath);
    end
end

nUnits = size(SingleUnits,2);
nStim = StimSettings.nStim;
nTypes = StimSettings.nTypes;
align2sniffs = StimSettings.align2sniffs;

if ~isempty(SniffsPlot)
    addSniffPlot = 1;
    AllSpikes(nUnits+1) = {SniffsPlot};
end

%% Actual Plotting
% general settings
switch StimSettings.SessionType
    case 'newCID'
        mycolors = brewermap(numel(nStim),'Dark2');
        trialsPerUnit = size(TTLs.Trial,1);
        plots_rows = 2;
        plots_cols = 6; %15;
        panelsPerPlot = plots_rows*plots_cols; % one per unit
        unitsPerPlot = panelsPerPlot;
        if bothPulses
            plotWidth = [-4 (StimSettings.timing(3)*5)/1000];
        else
            plotWidth = [-6 (StimSettings.timing(3)*2)/1000];
        end

    case '16_Odors'
        mycolors = brewermap(numel(nStim),'Dark2');
        trialsPerUnit = size(TTLs.Trial,1);
        plots_rows = 2;
        plots_cols = 5;
        panelsPerPlot = plots_rows*plots_cols; % one per unit
        unitsPerPlot = panelsPerPlot;
        plotWidth = (StimSettings.timing(3)*2)/1000*[-0.5 1];
        %plotWidth = [-6 (StimSettings.timing(3)*2)/1000];

    case '16_Concs*'
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        for c = 0:0.2:0.6
            lighter = basecolors + (1 - basecolors)*c;
            mycolors = vertcat(lighter, mycolors);
        end
        trialsPerUnit = size(TTLs.Trial,1);
        plots_rows = 4;
        plots_cols = 10;
        panelsPerPlot = plots_rows*plots_cols; % one per unit
        unitsPerPlot = panelsPerPlot;
        plotWidth = (StimSettings.timing(3)*2)/1000*[-0.5 1];
        %plotWidth = [-6 (StimSettings.timing(3)*2)/1000];
        StimSettings.SessionType = '16_Odors';

    case '16_Concs'
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        scale = 0.6:-0.2:0;
        for c = 1:4
            mycolors(:,:,c) = basecolors + (1 - basecolors)*scale(c);
        end
        mycolors = reshape(permute(mycolors, [3 1 2]), 20, 3);      % rows: [color1's N shades; color2's N shades; ...]
        trialsPerUnit = size(TTLs.Trial,1);
        plots_rows = 4;
        plots_cols = 10;
        panelsPerPlot = plots_rows*plots_cols; % one per unit
        unitsPerPlot = panelsPerPlot;
        plotWidth = (StimSettings.timing(3)*2)/1000*[-0.5 1];
        %plotWidth = [-6 (StimSettings.timing(3)*2)/1000];
        StimSettings.SessionType = '16_Odors';

    case 'otherwise'
        nConcs = numel(nTypes);
        extraColors = 2;
        mycolors = brewermap(2*(nConcs+extraColors),'YlOrRd');
        mycolors(1:(2*extraColors),:) = []; % too light
        trialsPerOdor = size(TTLs.Trial,1)/nStim;
        unitsPerPlot = 5;
        plots_rows = unitsPerPlot;
        plots_cols = numel(nStim);
        panelsPerPlot = plots_rows*plots_cols; % one per unit
        %plotWidth = [-6 (StimSettings.timing(3)*2)/1000];
        plotWidth = (StimSettings.timing(3)*2)/1000*[-0.5 1];
end
FigPosition = [2150 80 1750 900];

%% Plotting Sniffs
figure;

%% group by stimuli
subplot(1,4,1);
hold on;
odorON = [];
repsDone = 0;
for odor = 1:numel(nStim)
    for conc = 1:numel(nTypes) % this is just one
        SpikesPlot = [];
        whichtrials = intersect(find(TTLs.Trial(:,4)==nStim(odor)),find(TTLs.Trial(:,5)==nTypes(conc)));
        nreps = numel(whichtrials);
        for rep = 1:numel(whichtrials)
            thisTrial = whichtrials(rep);
            thistrialspikes = [];
            thistrialspikes = find(SniffsPlot(:,2)==thisTrial);
            thistrialspikes = SniffsPlot(thistrialspikes,1);
            SpikesPlot = vertcat(SpikesPlot, ...
                [thistrialspikes ...
                (rep + repsDone)*ones(numel(thistrialspikes),1)]...
                );
            if size(TTLs.Trial,2) == 12 && align2sniffs
                odorON = vertcat(odorON, [TTLs.Trial(thisTrial,12) , (rep + repsDone)]);
            else
                odorON = vertcat(odorON, [0 , (rep + repsDone)]);
            end
        end
        repsDone = repsDone + nreps;
        plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 4, 'color', mycolors(odor,:));
    end
end

plotHeight = (max(odorON(:,2))+1)/500;
set(gca,'XTick',[],'YTick',[],'XLim',plotWidth,'YLim',[0 plotHeight]);
plot(odorON(:,1),odorON(:,2)/500,'k');
plot(odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');
% plot second odor pulse if needed
if strcmp(StimSettings.SessionType,'newCID')
    pulseOffset = 1+ (StimSettings.timing(3) + StimSettings.timing(4))/1000;
    plot(pulseOffset+odorON(:,1),odorON(:,2)/500,'k');
    plot(pulseOffset+odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');
end

%% group by trial occurence, align to odor onset
TrialOrder = 1:size(TTLs.Trial,1);

subplot(1,4,2);
hold on;
[odorON] = PlotSniffsByTrialOrder(TrialOrder,TTLs,SniffsPlot,-1,plotWidth);
plot(odorON(:,1),odorON(:,2)/500,'r','LineWidth',0.5);

subplot(1,4,3);
hold on;
PlotSniffsByTrialOrder(TrialOrder,TTLs,SniffsPlot,0,plotWidth);

subplot(1,4,4);
hold on;
PlotSniffsByTrialOrder(TrialOrder,TTLs,SniffsPlot,2,plotWidth);

function [odorON] = PlotSniffsByTrialOrder(TrialOrder,TTLs,SniffsPlot,adjustcase,plotWidth)
odorON = [];
for myTrial = 1:numel(TrialOrder)
    thisTrial = TrialOrder(myTrial);
%     odor = find(nStim==TTLs.Trial(myTrial,4));
%     conc = find(nTypes==TTLs.Trial(myTrial,5));
    thistrialspikes = find(SniffsPlot(:,2)==thisTrial);
    thistrialspikes = SniffsPlot(thistrialspikes,1);
    switch adjustcase
        case 0
            adjustby = 0;
        case -1
            % align to odor onset
            adjustby = - TTLs.Trial(thisTrial,12);
        case 2
            % align to second sniff
            adjustby = -thistrialspikes(find(thistrialspikes>0,1,"first"));
        otherwise
            adjustby = 0;
    end
    SpikesPlot = [thistrialspikes myTrial*ones(numel(thistrialspikes),1)];
    odorON = vertcat(odorON, [TTLs.Trial(thisTrial,12)+adjustby , myTrial]);
    plot(SpikesPlot(:,1)+adjustby, SpikesPlot(:,2)/500, '.k','Markersize', 4, 'color', 'k');
    %plot(SpikesPlot(:,1)+adjustby, SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(odor,:));
end
plotHeight = (max(odorON(:,2))+1)/500;
set(gca,'XTick',[],'YTick',[],'XLim',plotWidth,'YLim',[0 plotHeight]);
end