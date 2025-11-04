%function [] = FiveOdorsSniffResponseExaminer() %(myDir,myStimFile)

%% function to plot sniff locked responses of individual neurons 
% aligned to sniffs and during odor presentatiion
% optimized for conc-identity experiments done by Marie

whichmouse = 'E6'; % 'E2', 'E6';
window = [-0.1 1];
binsize = 0.01;
savefigs = 1;
savePSTHs = 0;
unitsPerFig = 10;

switch whichmouse
    case 'E2' % AON
        sessionname = '2022-06-11_13-57-38_cid-processed.mat';
    case 'E3' % AON mice
        sessionname = '2022-06-14_11-36-14_cid-processed.mat';
    case 'E6' % APC mouse
        sessionname = '2022-06-10_11-40-39_cid-processed.mat';
end

FullProcessedPath = fullfile('/mnt/data/CID/Processed',whichmouse,sessionname);
FiguresPath = fullfile('/mnt/data/CID/Processed',whichmouse,'OdorMaps',strrep(sessionname,'cid-processed.mat','SniffOdorSummary.pdf'));

load(FullProcessedPath);

if savefigs
    if ~exist(fileparts(FiguresPath),'dir')
        mkdir(fileparts(FiguresPath));
    end
end

%% load sniffs and assign stimulus identity to each sniff
if exist('CuratedMFSSniffTimestamps','var')
    CuratedSniffTimestamps = CuratedMFSSniffTimestamps;
end
if ~exist('CuratedSniffTimestamps','var')
    keyboard;
else
% create AllSniffs
% cols that matter : 
% 1 (TS of inhalation start)
% 2 (Inhalation End)
% 3 (Inhalation Duration)
% 4 (Air on or off?)
% 5 (odor identity)
% 6 (odor conc)
% 12 (idx of Inhalatation start)
% 13 (idx of next sniff start)
% 
if ~any(abs(CuratedSniffTimestamps(1:end-1,3)-CuratedSniffTimestamps(2:end,1))>=0.003)
    AllSniffs(:,1:2) = CuratedSniffTimestamps(:,1:2);
    AllSniffs(:,3) = AllSniffs(:,2) - AllSniffs(:,1);
    AllSniffs(:,11:12) = CuratedSniffTimestamps(:,8:9);
    AllSniffs(1:end-1,13) = CuratedSniffTimestamps(2:end,8);
else
    keyboard;
end
end

% give each sniff an instantaneous freq - how close was the
% previous sniff (col 9), how close is the next sniff (col 10)
AllSniffs(:,9) = AllSniffs(:,1) - vertcat(nan, AllSniffs(1:end-1,1));
AllSniffs(:,10) = vertcat(AllSniffs(2:end,1), nan) - AllSniffs(:,1);

TTLs.Trial((TTLs.Trial(:,4)==0),4) =  max(TTLs.Trial(:,4)) + 1;
% add the column infos to the sniffs
for i = 1:size(TTLs.Trial,1) % every trial
    TS = TTLs.Trial(i,7:8); % odor1 start, stop, purge start, stop
    for j = 1:size(TS,1)
        whichsniffs = find( (AllSniffs(:,1)>=TS(j,1)) & (AllSniffs(:,1)<TS(j,2)) );
        AllSniffs(whichsniffs,5) = TTLs.Trial(i,4); % odor identity
        AllSniffs(whichsniffs,6) = TTLs.Trial(i,5);
    end
    % a hack assign sniffs to a given trial
    thisTrialsniffs = find( (AllSniffs(:,1)>=TTLs.Trial(i,1)) & (AllSniffs(:,1)<TTLs.Trial(i,2)) );
    AllSniffs(thisTrialsniffs,7) = TTLs.Trial(i,4);
end

StimSettings.Odors = unique(TTLs.Trial(:,4));
StimSettings.Concs = unique(TTLs.Trial(:,5));

%% group sniffs by odors
% there are air OFF sniffs here
AllSniffs(:,4) = 1; % air is always on
[ParsedSniffs, StimulusList] = ParseSniffsByStimuli(AllSniffs, 'SortBy', 5); % 5 - sort by session phase, location (conc) and then occurence

%% Load spikes
% SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

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
    
    % plot the 5 odors, start with highest concentration
    for conc = 1:numel(StimSettings.Concs)
        for odor = 1:numel(StimSettings.Odors)
            whichstim = find(StimulusList == StimSettings.Odors(odor));
            Sniffs2Use{1} = ParsedSniffs{whichstim+1};
            Sniffs2Use{1} = sortrows(Sniffs2Use{1},[6 1]);

            % select only sniffs that belong to the selected stim intensity
            Sniffs2Use{1} = Sniffs2Use{1}(find(Sniffs2Use{1}(:,6)==StimSettings.Concs(conc)),:);

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
                    FiguresPath,...
                    'ContentType','vector');
            else
                exportgraphics(gcf, ...
                    FiguresPath,...
                    'ContentType','vector','Append',true);
            end
        end
        close(gcf);
    end
end

if savePSTHs
    save(FullProcessedPath,"AllunitsPSTH","SniffChunks","SniffsUsed","durations",'-append');
end
%end