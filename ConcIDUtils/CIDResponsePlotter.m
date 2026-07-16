
%[StimSettings, TTLs, SingleUnits, AllSpikes, TrialWiseSniffs, SniffsPlot] = CIDResponsePrepper(myKsDir);

plots_rows = 3;
plots_cols = 12; %15;
%UnitList = [SingleUnits.id]';
%UnitList = [0 63 66 8 12 18 76 33 35 38 24 22 32 39 90 85 88 48 57 44 51 61 60 112 120 123 114 108 109 125 105 127];
UnitList =  [76 48 44 51 120 105 127 0 60 125 112 114 12 90 61 32 39 38 85 88 109 123 63 35 66 18 33 57 24 22 8 108];

% for saving figures to pdf
savefigs = 0;
if savefigs
    FigPath = strrep(myKsDir,'/mnt/data/Sorted','/home/priyanka/Desktop/cid');
    if ~exist("FigPath",'dir')
        mkdir(FigPath);
    end
end

nUnits = numel(UnitList);
nStim = StimSettings.nStim;
nTypes = StimSettings.nTypes;
align2sniffs = StimSettings.align2sniffs;

addSniffPlot = 0;
% if ~isempty(SniffsPlot)
%     addSniffPlot = 1;
%     AllSpikes(nUnits+1) = {SniffsPlot};
% end

%% Actual Plotting
% general settings
FigPosition = [2150 80 1750 900];
trialsPerUnit = size(TTLs.Trial,1);
panelsPerPlot = plots_rows*plots_cols; % one per unit
unitsPerPlot = panelsPerPlot;
if strcmp(StimSettings.SessionType,'newCID')
    if bothPulses
        plotWidth = [-4 (StimSettings.timing(3)*5)/1000];
    else
        plotWidth = [-6 (StimSettings.timing(3)*2)/1000];
    end
else
    plotWidth = (StimSettings.timing(3)*2)/1000*[-0.5 1];
end

switch StimSettings.SessionType
    case {'newCID', '16_Odors'}
        mycolors = brewermap(numel(nStim),'Dark2');

    case {'16_Concs*'}
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        for c = 0:0.2:0.6
            lighter = basecolors + (1 - basecolors)*c;
            mycolors = vertcat(lighter, mycolors);
        end
        StimSettings.SessionType = '16_Odors';

    case '16_Concs'
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        scale = 0.6:-0.2:0;
        for c = 1:4
            mycolors(:,:,c) = basecolors + (1 - basecolors)*scale(c);
        end
        mycolors = reshape(permute(mycolors, [3 1 2]), 20, 3);      % rows: [color1's N shades; color2's N shades; ...]
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

%% Actual Plotting
for m = 1:(nUnits+addSniffPlot)
    % Figure/subplot handling
    if mod(m,panelsPerPlot) == 1
        FigureName = ['Units ',num2str(m),'-',num2str(m+unitsPerPlot-1)];
        figure('Name',FigureName);
        set(gcf,'Position',FigPosition);
    end
    
    if m <= nUnits
        n = find([SingleUnits.id]==UnitList(m));
    else
        n = m;
    end
    % if 16 odors, loop by odors per unit for one subplot,
    if strcmp(StimSettings.SessionType,'newCID')|strcmp(StimSettings.SessionType,'16_Odors')
        whichplot = mod(m,panelsPerPlot);
        if ~whichplot
            whichplot = panelsPerPlot;
        end
        subplot(plots_rows,plots_cols,whichplot);
        hold on
        odorON = [];
        repsDone = 0;
        if m <= nUnits
            title(['Unit id: ',num2str(SingleUnits(n).id)]);
        else
            title('Sniff Raster');
        end
        for odor = 1:numel(nStim)
            for conc = 1:numel(nTypes) % this is just one
                SpikesPlot = [];
                whichtrials = intersect(find(TTLs.Trial(:,4)==nStim(odor)),find(TTLs.Trial(:,5)==nTypes(conc)));
                nreps = numel(whichtrials);
                for rep = 1:numel(whichtrials)
                    thisTrial = whichtrials(rep);
                    thistrialspikes = [];
                    thistrialspikes = find(AllSpikes{n}(:,2)==thisTrial);
                    thistrialspikes = AllSpikes{n}(thistrialspikes,1);
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
                plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(odor,:));
            end
        end

        plotHeight = (max(odorON(:,2))+1)/500;
        set(gca,'XTick',[],'YTick',[],'XLim',plotWidth,'YLim',[0 plotHeight]);
        plot(odorON(:,1),odorON(:,2)/500,'k');
        plot(odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');
        % plot second odor pulse if needed
        if strcmp(StimSettings.SessionType,'newCID')
            if bothPulses
                pulseOffset = 1+ (StimSettings.timing(3) + StimSettings.timing(4))/1000;
                plot(pulseOffset+odorON(:,1),odorON(:,2)/500,'k');
                plot(pulseOffset+odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');
            end
        end
        if m <= nUnits
            if SingleUnits(n).quality == 2
                set(gca,'Box','on');
            end
        end
    end

    % if conc. series, loop by conc. per odor per unit for one subplot
    if strcmp(StimSettings.SessionType,'ConcentrationSeries')
        whichRow = mod(m,unitsPerPlot);
        if ~whichRow
            whichRow = unitsPerPlot;
        end

        for odor = 1:numel(nStim)
            whichplot = ((whichRow-1)*plots_cols) + odor;
            subplot(plots_rows,plots_cols,whichplot);
            hold on
            odorON = [];
            repsDone = 0;

            for conc = 1:numel(nTypes) % this is just one
                SpikesPlot = [];
                whichtrials = intersect(find(TTLs.Trial(:,4)==nStim(odor)),find(TTLs.Trial(:,5)==nTypes(conc)));
                nreps = numel(whichtrials);
                for rep = 1:numel(whichtrials)
                    thisTrial = whichtrials(rep);
                    thistrialspikes = [];
                    thistrialspikes = find(AllSpikes{n}(:,2)==thisTrial);
                    thistrialspikes = AllSpikes{n}(thistrialspikes,1);
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
                plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(conc,:));
            end

            plotHeight = (max(odorON(:,2))+1)/500;
            set(gca,'XTick',[],'YTick',[],'XLim',plotWidth,'YLim',[0 plotHeight]);
            plot(odorON(:,1),odorON(:,2)/500,'k');
            plot(odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');

            if SingleUnits(n).quality == 2
                set(gca,'Box','on');
            end
        end
    end

    % if last panel for the figure
    if ~mod(m,panelsPerPlot)
        if savefigs
            saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
            close(gcf);
        end
    end
end
