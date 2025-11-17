function [] = SixteenOdorsSniffResponseExaminer() %(myDir,myStimFile)

whichmouse = 'T3';
window = [-0.1 1];
binsize = 0.01;
savefigs = 0;
unitsPerFig = 10;

switch whichmouse
    case 'T3'
        myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31/';
    case 'T2'
        myKsDir = '/mnt/data/Sorted/T2/_2025-05-21_09-18-56_2025-05-21_11-11-12_2025-05-21_11-23-51/';
end

%% load sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'),"SniffCoords");
if exist("SniffCoords")
    % make AllSniffs from SniffCoords
    switch whichmouse
        case 'T3'
            AllSniffs(:,1:2) = SniffCoords(:,1:2); % if using thermistor peaks
            AllSniffs(:,11:12) = SniffCoords(:,4:5); % if using thermistor peaks
        case 'T2'
            AllSniffs(:,1:2) = SniffCoords(:,6:7); % if using mfs2thermistor peaks
            AllSniffs(:,11:12) = SniffCoords(:,9:10); % if using mfs2thermistor peaks
    end
    AllSniffs(:,3)   = SniffCoords(:,2) - SniffCoords(:,1);
    AllSniffs(1:end-1,13) = AllSniffs(2:end,11); 
    AllSniffs(end,13) = AllSniffs(end,13) + 1; 

    % add odorTTL info
    load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings');
    
    % add the column infos to the sniffs
    for i = 1:size(TTLs.Trial,1) % every trial
        t = TTLs.Trial(i,7:10); % odor1 start, stop, purge start, stop
        t(end+1) = t(end) + TTLs.Trial(i,9)-TTLs.Trial(i,8); % add a post-purge period the same as the first pulse
        TS = [t(1:end-1)' t(2:end)' [1 3 2 4]'];
        for j = 1:size(TS,1)
            whichsniffs = find( (AllSniffs(:,1)>=TS(j,1)) & (AllSniffs(:,1)<TS(j,2)) );
            AllSniffs(whichsniffs,5) = TTLs.Trial(i,4); % odor identity
            AllSniffs(whichsniffs,6) = TS(j,3);
        end
        % a hack assign sniffs to a given trial
        thisTrialsniffs = find( (AllSniffs(:,1)>=TTLs.Trial(i,1)) & (AllSniffs(:,1)<TTLs.Trial(i,2)) );
        AllSniffs(thisTrialsniffs,7) = TTLs.Trial(i,4);
    end
else
    %%
    if exist(fullfile(myKsDir,'quickprocesssniffs.mat'))
        load(fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs', 'RespirationData');
    end
end

%% hack 2: 
% give each sniff an instantaneous freq - how close was the
% previous sniff (col 9), how close is the next sniff (col 10)
AllSniffs(:,9) = AllSniffs(:,1) - vertcat(nan, AllSniffs(1:end-1,1));
AllSniffs(:,10) = vertcat(AllSniffs(2:end,1), nan) - AllSniffs(:,1);

%% separate session phases if necessary
load(fullfile(myKsDir,'SessionDetails.mat'));
startTS = 0;
for i = 1:numel(Files.Samples)
    endTS = Files.Samples(i)/30000 + startTS; % in seconds
    whichsniffs = find( (AllSniffs(:,1)>=startTS) & (AllSniffs(:,1)<endTS) );
    AllSniffs(whichsniffs,8) = i; % session phase
    startTS = endTS;
end

%% group sniffs by odors
% there are air OFF sniffs here
AllSniffs(:,4) = 1; % air is always on
[ParsedSniffs, StimulusList] = ParseSniffsByStimuli(AllSniffs, 'SortBy', 5); % 5 - sort by session phase, location (conc) and then occurence

%% Load spikes
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

%% load odor stimulus info
if exist(fullfile(myKsDir,'quickprocessOdorTTLs.mat'))
    load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'StimSettings');
end

%% figures related

nRows = 1; nCols = unitsPerFig + 1;
mycolors = brewermap(18,'Accent');
durations = [];
AllunitsPSTH = [];
SniffsUsed = [];
SniffChunks = [];

for n = 1:nUnits
    if rem(n,unitsPerFig) == 1
        FigureName = ['Unit ',num2str(n),' to Unit ',num2str(n+unitsPerFig-1)]; % one figure per cell
        figure('Name',FigureName);
        set(gcf,'Position',[2000   138   1366   700]);
    end
    if rem(n,unitsPerFig) == 0
        whichsubplot = unitsPerFig;
    else
        whichsubplot = rem(n,unitsPerFig);
    end
    subplot(nRows,nCols,whichsubplot);
    hold on;
    sniffsDone = 0;
    mypsth = [];

    thisUnitSpikes = SingleUnits(n).spikes;
    
    % plot the air spikes
    Sniffs2Use{1} = ParsedSniffs{2}; % Air sniffs
    Sniffs2Use{1}(:,8) = 1;

    % select only air sniffs that occured before the first trial
    untilSniff = find(Sniffs2Use{1}(:,1)>=TTLs.Trial(1,1),1,'first');
    Sniffs2Use{1}(untilSniff:end,:) = [];

    if n == 1
        durations = Sniffs2Use{1}(:,9:10);
        sniffgrouping = round(linspace(1,size(Sniffs2Use{1},1),5));
        sniffgrouping(2,1:end-1) = sniffgrouping(1,2:end) - 1;
        sniffgrouping = sniffgrouping';
        sniffgrouping(end,:) = [];
        SniffsUsed = vertcat(SniffsUsed, Sniffs2Use{1});
        SniffChunks = vertcat(SniffChunks, size(Sniffs2Use{1},1));
    end

    % get sniff aligned spikes - default window is -0.1 to 1
    %[SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    [SpikeRaster, maxsniffs, PSTHout] = GetSniffLockedSpikesAndPSTH(Sniffs2Use, thisUnitSpikes, 'window', window, 'binsize', binsize);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    sniffsDone = sniffsDone + maxsniffs;

    % get airPSTHs
    for q = 1:size(sniffgrouping,1)
        chunkedsniffs{1} = Sniffs2Use{1}(sniffgrouping(q,1):sniffgrouping(q,2),:);
        [~, ~, PSTHout] = GetSniffLockedSpikesAndPSTH(chunkedsniffs, thisUnitSpikes, 'window', window, 'binsize', binsize);
        mypsth = vertcat(mypsth, PSTHout{1});
    end
    
    % plot the 10 mini odors
    for odor = 1:numel(StimSettings.miniOdors)
        whichstim = find(StimulusList == StimSettings.miniOdors(odor));
        Sniffs2Use{1} = ParsedSniffs{whichstim+1};
        Sniffs2Use{1} = sortrows(Sniffs2Use{1},[6 1]);
        
        % select only sniffs that occured in first stimulus (not the purge)
        untilSniff = find(Sniffs2Use{1}(:,6)>2,1,'first');
        Sniffs2Use{1}(untilSniff:end,:) = [];

        %[SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes, sniffsDone);
        [SpikeRaster, maxsniffs, PSTHout] = GetSniffLockedSpikesAndPSTH(Sniffs2Use, thisUnitSpikes, 'yoffset', sniffsDone, 'window', window, 'binsize', binsize);
        mypsth = vertcat(mypsth, PSTHout{1});
        plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
        
        % odor Line
        plot([-0.15 -0.15],sniffsDone+[1 maxsniffs],'LineWidth',4,'Color',mycolors(whichstim+1,:));
        sniffsDone = sniffsDone + maxsniffs;
        
        if n == 1
            durations = vertcat(durations, Sniffs2Use{1}(:,9:10));
            SniffsUsed = vertcat(SniffsUsed, Sniffs2Use{1});
            SniffChunks = vertcat(SniffChunks, size(Sniffs2Use{1},1));
        end
    end
    
    % plot the 5 mega odors - skip the blank
    for odor = 1:numel(StimSettings.megaOdors)
        whichstim = find(StimulusList == StimSettings.megaOdors(odor));
        Sniffs2Use{1} = ParsedSniffs{whichstim+1};
        Sniffs2Use{1} = sortrows(Sniffs2Use{1},[6 1]);

        % select only sniffs that occured in first stimulus (not the purge)
        untilSniff = find(Sniffs2Use{1}(:,6)>2,1,'first');
        Sniffs2Use{1}(untilSniff:end,:) = [];

        %[SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes, sniffsDone);
        [SpikeRaster, maxsniffs, PSTHout] = GetSniffLockedSpikesAndPSTH(Sniffs2Use, thisUnitSpikes, 'yoffset', sniffsDone, 'window', window, 'binsize', binsize);
        plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
        mypsth = vertcat(mypsth, PSTHout{1});

        % odor Line
        plot([-0.15 -0.15],sniffsDone+[1 maxsniffs],'LineWidth',4,'Color',mycolors(whichstim+1,:));
        sniffsDone = sniffsDone + maxsniffs;
        
        if n == 1
            durations = vertcat(durations, Sniffs2Use{1}(:,9:10));
            SniffsUsed = vertcat(SniffsUsed, Sniffs2Use{1});
            SniffChunks = vertcat(SniffChunks, size(Sniffs2Use{1},1));
        end
    end
    
    AllunitsPSTH{n} = mypsth;

    if n == 1
        b = get(gca,'YLim');
        durations(end+1:b(2),:) = nan;
    end

    set(gca,'TickDir','out','YTick',[]);

    if rem(n,unitsPerFig)==0 || n == nUnits
        subplot(nRows,nCols,whichsubplot+1);
        imagesc(flipud(durations),[0.1 0.5])
        colormap(brewermap([100],'PiYg'));
        colorbar;
        set(gca,'TickDir','out','YTick',[]);
    end

    if (rem(n,unitsPerFig)==0 || n == nUnits)
        if savefigs
            set(gcf,'Color','w');
            set(gcf,'renderer','Painters');
            if n == unitsPerFig
                exportgraphics(gcf, ...
                    fullfile(myKsDir,'OdorMaps','SniffOdorSummaryFirstPulse.pdf'),...
                    'ContentType','vector');
            else
                exportgraphics(gcf, ...
                    fullfile(myKsDir,'OdorMaps','SniffOdorSummaryFirstPulse.pdf'),...
                    'ContentType','vector','Append',true);
            end
        end
        close(gcf);
    end
end

save(fullfile(myKsDir,'OdorMaps','SniffOdorSummaryFirstPulseParams.mat'),"AllunitsPSTH","SniffChunks","SniffsUsed","durations");

end