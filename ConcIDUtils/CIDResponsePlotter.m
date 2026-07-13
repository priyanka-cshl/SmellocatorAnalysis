%% Input
%myKsDir = '/mnt/storage/Sorted/Q8/2022-12-19_15-41-55';
%myKsDir = '/media/priyanka/ABC-ntfs/EphysSorted/Q9/2022-12-15_16-28-22';

%% Plotting/Figure related settings
savefigs = 0;
align2sniffs = 0; % first sniff, second sniff
whichSniffSensor = 2; % 1 = thermistor, 2 = MFS
stackUnits = 0;
stackConcs = 1;
% plotStyle = 'sniffwise';
if savefigs
    FigPath = strrep(myKsDir,'/mnt/data/Sorted','/home/priyanka/Desktop/cid');
    if ~exist("FigPath",'dir')
        mkdir(FigPath);
    end
end

%% load the data
clear KS4Units;
load(fullfile(myKsDir,'quickprocesssniffs.mat')); % sniff times, KS4Units
load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings'); % odorTTLs
SingleUnits = LoadKS4Units(myKsDir,'minSpikes',0.25);
disp(['found ',num2str(size(SingleUnits,2)),' good, >0.25Hz units']);
nUnits = size(SingleUnits,2);

%% stimulus settings for plotting
if stackConcs
    if strcmp(StimSettings.SessionType,'ConcentrationSeries')
        newStimCol = (TTLs.Trial(:,5)*10^4) + TTLs.Trial(:,4);
        newTypesCol = TTLs.Trial(:,5)*0;
        TTLs.Trial(:,4:5) = [newStimCol newTypesCol];
        StimSettings.SessionType = '16_Concs';
    end
end
nStim = unique(TTLs.Trial(:,4)); % no. of odors delivered
nTypes = unique(TTLs.Trial(:,5)); % concentrations used
%nreps = max(TTLs.Trial(:,6));

%% Hacks for some experiments
% make air trials as last stimulus
TTLs.Trial((TTLs.Trial(:,4)==0),4) =  max(TTLs.Trial(:,4)) + 1;

if any(nStim == -1)
    f = find(TTLs.Trial(:,4)==-1);
    TTLs.Trial(f,5) = nTypes(2);
    nTypes(1,:) = [];
    nStim(nStim<0,:) = [];
    TTLs.Trial(f,7) = TTLs.Trial(f,1) + StimSettings.timing(2)/1000;
    TTLs.Trial(f,8) = TTLs.Trial(f,7) + StimSettings.timing(3)/1000;
end

% for new CID experiments
if size(TTLs.Odor,1) == 2*size(TTLs.Trial,1)
    StimSettings.SessionType = 'newCID';
end

if ~isfield(StimSettings,'SessionType')
    keyboard;
end

%% process sniffs
if align2sniffs
    % find first inhalation start after odor start
    if whichSniffSensor==1 && exist('CuratedSniffTimestamps','var')
        MySniffTimeStamps = CuratedSniffTimestamps(:,1:3);
    elseif whichSniffSensor==2 && exist('CuratedMFSSniffTimestamps','var')
        MySniffTimeStamps = CuratedMFSSniffTimestamps(:,1:3);
    else
        keyboard;
    end
    TrialWiseSniffs = [];
    for t = 1:size(TTLs.Trial,1)
        if TTLs.Trial(t,4)>0 % every valid trial
            ts = TTLs.Trial(t,[1 2 7 8]); % trial start, stop, odor start, stop
            % find all sniffs
            thisTrialSniffs = intersect(find(MySniffTimeStamps(:,1)>=ts(1)),find(MySniffTimeStamps(:,1)<ts(2)));
            firstsniff = find(MySniffTimeStamps(thisTrialSniffs,1)>=ts(3),1,'first');
            TTLs.Trial(t,10) = MySniffTimeStamps(thisTrialSniffs(firstsniff+(align2sniffs-1)),1); % time of first sniff (true time)
            TTLs.Trial(t,11) = TTLs.Trial(t,7) - TTLs.Trial(t,10); % % time of first sniff from odor Onset
            thisTrialSniffs = MySniffTimeStamps(thisTrialSniffs,1:3) - ts(3); % odor start
            % now lets index every sniff
            thisTrialSniffs(:,4) = t;
            aftersniffs = find(thisTrialSniffs(:,1)>=0);
            beforesniffs = 1:(aftersniffs(1)-1);
            thisTrialSniffs(beforesniffs,5) = -flipud(beforesniffs(:));
            thisTrialSniffs(aftersniffs,5)  = (1:numel(aftersniffs))-1;
            thisTrialSniffs(:,6) = ts(3); % useful to get back actual value
            TrialWiseSniffs = vertcat(TrialWiseSniffs, thisTrialSniffs);
        end
    end
end

% TTLs.Trial has now sniff info in additional columns
% 10    : true time of first (or nth) sniff after odor ON
% 11    : relative time of first (or nth) sniff after odor ON

% TrialWiseSniffs: Columns
% 1 to 3: inh start, end, next w.r.t. this Trial's odor Onset
% 4     : trial index
% 5     : sniff index within a trial, 0 = first sniff after odor onset
% 6     : actual odor ON time

%% Make  trial-aligned Spike Plot
for n = 1:nUnits
    SpikesPlot = [];
    thisUnitSpikes = SingleUnits(n).spikes;
    for t = 1:size(TTLs.Trial,1)
        if TTLs.Trial(t,4)>0
            % every trial
            ts = TTLs.Trial(t,[1 2 7 8]); % trial start, stop, odor start, stop
            if size(TTLs.Trial,2) >= 10 && align2sniffs
                ts(5) = TTLs.Trial(t, 10);
            else
                ts(5) = ts(3);
            end
            thistrialspikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<ts(2)));
            thistrialspikes = thisUnitSpikes(thistrialspikes) - ts(5);
            multiplier = t;
            SpikesPlot = vertcat(SpikesPlot, ...
                [thistrialspikes multiplier*ones(numel(thistrialspikes),1)]);
        end
    end
    AllSpikes(n) = {SpikesPlot};
end

%% Actual Plotting
if ~stackUnits
    % general settings
    switch StimSettings.SessionType
        case 'newCID'
            mycolors = brewermap(numel(nStim),'Dark2');
            trialsPerUnit = size(TTLs.Trial,1);
            plots_rows = 2;
            plots_cols = 6; %15;
            panelsPerPlot = plots_rows*plots_cols; % one per unit
            unitsPerPlot = panelsPerPlot;
            % plotWidth = [-6 (StimSettings.timing(3)*4)/1000];
            plotWidth = [-6 (StimSettings.timing(3)*2)/1000];

        case '16_Odors'
            mycolors = brewermap(numel(nStim),'Dark2');
            trialsPerUnit = size(TTLs.Trial,1);
            plots_rows = 2;
            plots_cols = 5;
            panelsPerPlot = plots_rows*plots_cols; % one per unit
            unitsPerPlot = panelsPerPlot;
            plotWidth = (StimSettings.timing(3)*2)/1000*[-0.5 1];
            %plotWidth = [-6 (StimSettings.timing(3)*2)/1000];

        case '16_Concs'
            mycolors = brewermap(5,'Dark2');
            for c = 1:4
                mycolors = vertcat(mycolors + (1 - mycolors)*c*0.2, ...
                    mycolors);
            end
            trialsPerUnit = size(TTLs.Trial,1);
            plots_rows = 4;
            plots_cols = 5;
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

    for n = 1:nUnits
        % Figure/subplot handling
        if mod(n,panelsPerPlot) == 1
            FigureName = ['Units ',num2str(n),'-',num2str(n+unitsPerPlot-1)];
            figure('Name',FigureName);
            set(gcf,'Position',FigPosition);
        end

        % if 16 odors, loop by odors per unit for one subplot,
        if strcmp(StimSettings.SessionType,'newCID')|strcmp(StimSettings.SessionType,'16_Odors')
            whichplot = mod(n,panelsPerPlot);
            if ~whichplot
                whichplot = panelsPerPlot;
            end
            subplot(plots_rows,plots_cols,whichplot);
            hold on
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
                        thistrialspikes = find(AllSpikes{n}(:,2)==thisTrial);
                        thistrialspikes = AllSpikes{n}(thistrialspikes,1);
                        SpikesPlot = vertcat(SpikesPlot, ...
                            [thistrialspikes ...
                            (rep + repsDone)*ones(numel(thistrialspikes),1)]...
                            );
                        if size(TTLs.Trial,2) == 11 && align2sniffs
                            odorON = vertcat(odorON, [TTLs.Trial(thisTrial,11) , (rep + repsDone)]);
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
        end

        % if conc. series, loop by conc. per odor per unit for one subplot
        if strcmp(StimSettings.SessionType,'ConcentrationSeries')
            whichRow = mod(n,unitsPerPlot);
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
                        if size(TTLs.Trial,2) == 11 && align2sniffs
                            odorON = vertcat(odorON, [TTLs.Trial(thisTrial,11) , (rep + repsDone)]);
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
            end
        end

        % if last panel for the figure
        if ~mod(n,panelsPerPlot)
            if savefigs
                saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                close(gcf);
            end
        end
    end

else % stackUnits
end
