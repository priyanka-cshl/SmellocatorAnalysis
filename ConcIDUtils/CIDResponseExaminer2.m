% general settings
savefigs = 0;
align2sniffs = 1; % first sniff, second sniff
stackUnits = 0;
plotStyle = 'sniffwise';

% load the data
clear KS4Units;
%myKsDir = '/mnt/storage/Sorted/Q8/2022-12-19_15-41-55';
%myKsDir = '/media/priyanka/ABC-ntfs/EphysSorted/Q9/2022-12-15_16-28-22';
load(fullfile(myKsDir,'quickprocesssniffs.mat')); % sniff times, KS4Units
load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings'); % odorTTLs

if exist('KS4Units') && isempty(dir(fullfile(myKsDir,'kilosort4','cluster_info*')))
    myUnits = [[KS4Units.id]' [KS4Units.tetrode]' [KS4Units.quality]'];
    myUnits(:,4) = 1:size(myUnits,1);
    % session wasn't curated in phy, keep only 'good' units
    myUnits(find(myUnits(:,3)~=2),:) = [];
    disp(['found ',num2str(size(myUnits,1)),' good units']);
    SingleUnits = KS4Units(myUnits(:,4));
else
    KS4Units = GetSingleUnits(fullfile(myKsDir, 'kilosort4'));
    myUnits = [[KS4Units.id]' [KS4Units.tetrode]' [KS4Units.quality]'];
    myUnits(:,4) = 1:size(myUnits,1);
    % session wasn't curated in phy, keep only 'good' units
    myUnits(find(myUnits(:,3)~=2),:) = [];
    disp(['found ',num2str(size(myUnits,1)),' good units']);
    SingleUnits = KS4Units(myUnits(:,4));
end

TTLs.Trial((TTLs.Trial(:,4)==0),4) =  max(TTLs.Trial(:,4)) + 1;

% FigPath = '/home/priyanka/Desktop/cid/E2/20220611';

%% stimulus settings for plotting
nStim = unique(TTLs.Trial(:,4));
nTypes = unique(TTLs.Trial(:,5));
nreps = max(TTLs.Trial(:,6));
nUnits = size(SingleUnits,2);

%StimSettings.SessionType = 'AllTrials';
if any(nStim == -1)
    f = find(TTLs.Trial(:,4)==-1);
    TTLs.Trial(f,5) = nTypes(2);
    nTypes(1,:) = [];
    TTLs.Trial(f,7) = TTLs.Trial(f,1) + StimSettings.timing(2)/1000;
    TTLs.Trial(f,8) = TTLs.Trial(f,7) + StimSettings.timing(3)/1000;
end

%% process sniffs
if align2sniffs
    % find first inhalation start after odor start
    if exist('CuratedSniffTimestamps','var')
        MySniffTimeStamps = CuratedSniffTimestamps(:,1:3);
    elseif exist('CuratedMFSSniffTimestamps','var')
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
            firstsniff = find(MySniffTimeStamps(thisTrialSniffs,1)>=TTLs.Trial(t,7),1,'first');
            TTLs.Trial(t,10) = MySniffTimeStamps(thisTrialSniffs(firstsniff+align2sniffs),1);
            TTLs.Trial(t,11) = TTLs.Trial(t,7) - TTLs.Trial(t,10);
            thisTrialSniffs = MySniffTimeStamps(thisTrialSniffs,1:3) - ts(3); % odor start
            % now lets index every sniff
            thisTrialSniffs(:,4) = t;
            aftersniffs = find(thisTrialSniffs(:,1)>=0);
            beforesniffs = 1:(aftersniffs(1)-1);
            thisTrialSniffs(beforesniffs,5) = -flipud(beforesniffs(:));
            thisTrialSniffs(aftersniffs,5)  = (1:numel(aftersniffs))-1;
            TrialWiseSniffs = vertcat(TrialWiseSniffs, thisTrialSniffs);
        end
    end
end

%%
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

%%

switch StimSettings.SessionType
    case 'ConcentrationSeries'
        %%
        if stackUnits
            mycolors = brewermap(16,'Dark2');
            nTrials = numel(StimSettings.Odors)*StimSettings.Reps;
            StimTable = TTLs.Trial(1:nTrials,4:6);

            % one plot per odor
            % one row per conc
            % one col per repeat (or stack repeats)

            odors = unique(StimTable(:,1));
            concs = unique(StimTable(:,2));
            reps  = unique(StimTable(:,3));

            for i = 1:numel(odors)
                figure;
                whichodor = odors(i);
                for j = 1:numel(concs)
                    whichconc = concs(j);
                    whichtrials = intersect(find(StimTable(:,1)==whichodor),find(StimTable(:,2)==whichconc));
                    % stack all repeats
                    RepOffset = linspace(0,1,1+numel(whichtrials));

                    SpikesPlot = [];
                    for n = 1:nUnits
                        for k = 1:numel(whichtrials)
                            thisTrial = whichtrials(k);
                            thistrialspikes = [];
                            thistrialspikes = find(AllSpikes{n}(:,2)==thisTrial);
                            thistrialspikes = AllSpikes{n}(thistrialspikes,1);
                            SpikesPlot = vertcat(SpikesPlot, ...
                                [thistrialspikes ...
                                ((n-1)+RepOffset(k))*ones(numel(thistrialspikes),1)]...
                                );
                        end
                    end
                    subplot(1,4,j);
                    plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2,'color', mycolors(i,:));
                    set(gca,'XTick',[],'YTick',[],'XLim',[-4 4])
                end
            end
            %%
        else
            mycolors = brewermap(2*(numel(find(nTypes~=0))+1),'YlOrRd');
            mycolors(1:2,:) = [];

            for n = 1:nUnits
                if mod(n,10) == 1
                    FigureName = ['Units ',num2str(n),'-',num2str(n+9)];
                    figure('Name',FigureName);
                end
                row = mod(n,10);
                if ~row
                    row = 10;
                end

                for odor = 1:numel(nStim)
                    if nStim(odor)>0
                        subplot(10,numel(nStim),((row-1)*numel(nStim) + odor));
                        hold on
                        odorON = [];
                        for conc = 1:numel(nTypes)
                            SpikesPlot = [];
                            whichtrials = intersect(find(TTLs.Trial(:,4)==nStim(odor)),find(TTLs.Trial(:,5)==nTypes(conc)));
                            for rep = 1:numel(whichtrials)
                                thisTrial = whichtrials(rep);
                                thistrialspikes = [];
                                thistrialspikes = find(AllSpikes{n}(:,2)==thisTrial);
                                thistrialspikes = AllSpikes{n}(thistrialspikes,1);
                                SpikesPlot = vertcat(SpikesPlot, ...
                                    [thistrialspikes ...
                                    (rep + (conc-1)*nreps)*ones(numel(thistrialspikes),1)]...
                                    );
                                if size(TTLs.Trial,2) == 11 && align2sniffs
                                    odorON = vertcat(odorON, [TTLs.Trial(thisTrial,11) , (rep + (conc-1)*nreps)]);
                                else
                                    odorON = vertcat(odorON, [0 , (rep + (conc-1)*nreps)]);
                                end
                            end
                            plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(conc*2,:));
                        end
                        set(gca,'XTick',[],'YTick',[],'XLim',[-6 4],'YLim',[0 0.056]);
                        plot(odorON(:,1),odorON(:,2)/500,'k');
                        plot(odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');
                    end
                end

                if ~mod(n,10)
                    set(gcf,'Position',[2150 80 1350 900]);
                    if savefigs
                        saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                        close(gcf);
                    end
                end

            end

            % last figure
            if mod(n,10) ~= 0
                set(gcf,'Position',[2150 80 1350 900]);
                if savefigs
                    saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                    close(gcf);
                end
            end

        end
        %%
    case '16_Odors'
        %%
        mycolors = brewermap(16,'Dark2');
        switch plotStyle
            case 'stackUnits'

            case 'sniffwise'
                for n = 1:nUnits
                    if mod(n,10) == 1
                        FigureName = ['Units ',num2str(n),'-',num2str(n+9)];
                        figure('Name',FigureName);
                    end
                    whichplot = mod(n,10);
                    if ~whichplot
                        whichplot = 10;
                    end

                    subplot(2,5,whichplot);
                    hold on
                    odorON = [];

                    for odor = 1:numel(nStim)
                        for conc = 1:numel(nTypes)
                            SpikesPlot = [];
                            whichtrials = intersect(find(TTLs.Trial(:,4)==nStim(odor)),find(TTLs.Trial(:,5)==nTypes(conc)));
                            for rep = 1:numel(whichtrials)
                                thisTrial = whichtrials(rep);
                                thistrialspikes = [];
                                thistrialspikes = find(AllSpikes{n}(:,2)==thisTrial);
                                thistrialspikes = AllSpikes{n}(thistrialspikes,1);

                                thistrialsniffs = 

                                SpikesPlot = vertcat(SpikesPlot, ...
                                    [thistrialspikes ...
                                    (rep + (odor-1)*nreps)*ones(numel(thistrialspikes),1)]...
                                    );
                                if size(TTLs.Trial,2) == 11 && align2sniffs
                                    odorON = vertcat(odorON, [TTLs.Trial(thisTrial,11) , (rep + (odor-1)*nreps)]);
                                else
                                    odorON = vertcat(odorON, [0 , (rep + (odor-1)*nreps)]);
                                end
                            end
                            plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(odor,:));
                        end
                    end

                    set(gca,'XTick',[],'YTick',[],'XLim',[-6 4],'YLim',[0 0.226]);
                    plot(odorON(:,1),odorON(:,2)/500,'k');
                    plot(odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');

                    if ~mod(n,10)
                        set(gcf,'Position',[2150 80 1350 900]);
                        if savefigs
                            saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                            close(gcf);
                        end
                    end
                end

                % last figure
                if mod(n,10) ~= 0
                    set(gcf,'Position',[2150 80 1350 900]);
                    if savefigs
                        saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                        close(gcf);
                    end
                end

            otherwise
                for n = 1:nUnits
                    if mod(n,10) == 1
                        FigureName = ['Units ',num2str(n),'-',num2str(n+9)];
                        figure('Name',FigureName);
                    end
                    whichplot = mod(n,10);
                    if ~whichplot
                        whichplot = 10;
                    end

                    subplot(2,5,whichplot);
                    hold on
                    odorON = [];

                    for odor = 1:numel(nStim)
                        for conc = 1:numel(nTypes)
                            SpikesPlot = [];
                            whichtrials = intersect(find(TTLs.Trial(:,4)==nStim(odor)),find(TTLs.Trial(:,5)==nTypes(conc)));
                            for rep = 1:numel(whichtrials)
                                thisTrial = whichtrials(rep);
                                thistrialspikes = [];
                                thistrialspikes = find(AllSpikes{n}(:,2)==thisTrial);
                                thistrialspikes = AllSpikes{n}(thistrialspikes,1);
                                SpikesPlot = vertcat(SpikesPlot, ...
                                    [thistrialspikes ...
                                    (rep + (odor-1)*nreps)*ones(numel(thistrialspikes),1)]...
                                    );
                                if size(TTLs.Trial,2) == 11 && align2sniffs
                                    odorON = vertcat(odorON, [TTLs.Trial(thisTrial,11) , (rep + (odor-1)*nreps)]);
                                else
                                    odorON = vertcat(odorON, [0 , (rep + (odor-1)*nreps)]);
                                end
                            end
                            plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(odor,:));
                        end
                    end

                    set(gca,'XTick',[],'YTick',[],'XLim',[-6 4],'YLim',[0 0.226]);
                    plot(odorON(:,1),odorON(:,2)/500,'k');
                    plot(odorON(:,1)+StimSettings.timing(3)/1000,odorON(:,2)/500,'k');

                    if ~mod(n,10)
                        set(gcf,'Position',[2150 80 1350 900]);
                        if savefigs
                            saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                            close(gcf);
                        end
                    end
                end

                % last figure
                if mod(n,10) ~= 0
                    set(gcf,'Position',[2150 80 1350 900]);
                    if savefigs
                        saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                        close(gcf);
                    end
                end
                %%


        end

    case 'ByRepeats'
        if numel(nStim) == 6
            maps = {'Greys','Blues','Oranges','Greens','Purples','Reds'};
        else
            maps = {'Blues','Oranges','Greens','Purples','Reds'};
        end

        for c = 1:numel(nStim)
            colors = brewermap(6,maps{c});
            mycolors(c) = {colors};
        end

        ureps = flipud(unique(TTLs.Trial(:,6)));
        trialsperrep = 20;
        for n = 1:nUnits
            if mod(n,10) == 1
                FigureName = ['Units ',num2str(n),'-',num2str(n+9)];
                figure('Name',FigureName);
            end
            whichplot = mod(n,10);
            if ~whichplot
                whichplot = 10;
            end
            thisUnitSpikes = SingleUnits(n).spikes;
            subplot(2,5,whichplot);
            hold on
            odorON = [];
            for thisrep = 1:numel(ureps)
                SpikesPlot = [];
                whichtrials = find(TTLs.Trial(:,6)==ureps(thisrep));
                for rep = 1:numel(whichtrials)
                    ts = TTLs.Trial(whichtrials(rep),[1 2 7 8]); % trial start, stop, odor start, stop
                    if size(TTLs.Trial,2) == 10
                        ts(5) = TTLs.Trial(whichtrials(rep), 10);
                    else
                        ts(5) = ts(3);
                    end
                    thistrialspikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<ts(2)));
                    thistrialspikes = thisUnitSpikes(thistrialspikes) - ts(5);
                    odorON = vertcat(odorON, [ts(3)-ts(5) , (rep + (thisrep-1)*trialsperrep)]);
                    SpikesPlot = vertcat(SpikesPlot, [thistrialspikes (rep + (thisrep-1)*trialsperrep)*ones(numel(thistrialspikes),1)]); % conc*ones(numel(thistrialspikes),1)]);
                end
                plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors{thisrep}(end,:));
            end
            set(gca,'XTick',[],'YTick',[],'XLim',[-6 4],'YLim',[0 0.4]);
            %         line([0 0],[0 0.056],'color','k');
            %         line([2 2],[0 0.056],'color','k');
            plot(odorON(:,1),odorON(:,2)/500,'k');
            plot(odorON(:,1)+2,odorON(:,2)/500,'k');

            if ~mod(n,10)
                set(gcf,'Position',[2150 80 1350 900]);
                if savefigs
                    saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                    close(gcf);
                end
            end
        end

        % last figure
        if mod(n,10) ~= 0
            set(gcf,'Position',[2150 80 1350 900]);
            if savefigs
                saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
                close(gcf);
            end
        end


    case 'AllOdors'





end