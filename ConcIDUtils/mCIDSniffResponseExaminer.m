function [] = mCIDSniffResponseExaminer() %(myDir,myStimFile)

myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31/';

%% load sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'),"SniffCoords");
if exist("SniffCoords")
    % make AllSniffs from SniffCoords
    AllSniffs(:,1:2) = SniffCoords(:,1:2); % if using thermistor peaks
    AllSniffs(:,11:12) = SniffCoords(:,4:5); % if using thermistor peaks
    AllSniffs(:,3)   = SniffCoords(:,2) - SniffCoords(:,1);
    AllSniffs(1:end-1,13) = AllSniffs(2:end,11); 
    AllSniffs(end,13) = AllSniffs(end,13) + 1; 

    % add odorTTL info
    load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings');
    
    % add the column infos to the sniffs
    for i = 1:size(TTLs.Trial,1) % everty trial
        t = TTLs.Trial(i,7:10); % odor1 start, stop, purge start, stop
        t(end+1) = t(end) + TTLs.Trial(i,9)-TTLs.Trial(i,8); % add a post-purge period the same as the first pulse
        TS = [t(1:end-1)' t(2:end)' [1 3 2 4]'];
        for j = 1:size(TS,1)
            whichsniffs = find( (AllSniffs(:,1)>=TS(j,1)) & (AllSniffs(:,1)<TS(j,2)) );
            AllSniffs(whichsniffs,5) = TTLs.Trial(i,4); % odor identity
            AllSniffs(whichsniffs,6) = TS(j,3);
        end
    end
else
    %%
    if exist(fullfile(myKsDir,'quickprocesssniffs.mat'))
        load(fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs', 'RespirationData');
    end
end

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
[ParsedSniffs, StimulusList] = ParseSniffsByStimuli(AllSniffs, 'SortBy', 3);

%% Load spikes
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

%% load odor stimulus info
if exist(fullfile(myKsDir,'quickprocessOdorTTLs.mat'))
    load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'StimSettings');
end

%% figures related
savefigs = 0;
% mycolors = brewermap(10,'YlOrRd');
% mycolors(1:2,:) = [];
nRows = 4; nCols = 6;

for n = 1:nUnits
    FigureName = ['Unit ',num2str(n)]; % one figure per cell
    figure('Name',FigureName);
    thisUnitSpikes = SingleUnits(n).spikes;
    
    % plot the air spikes
    Sniffs2Use{1} = ParsedSniffs{2}; % Air sniffs
    Sniffs2Use{1}(:,8) = 1;
    % resort if necessary
    
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    whichsubplots = 1:nCols:(nRows*nCols);
    subplot(nRows,nCols,whichsubplots);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    set(gca,'XTick',[],'YTick',[]);
    
    % plot the 10 mini odors
    for odor = 1:numel(StimSettings.miniOdors)
        if odor > 5
            subplot(nRows,nCols,odor+2);
        else
            subplot(nRows,nCols,odor+1);
        end
        
        whichstim = find(StimulusList == StimSettings.miniOdors(odor));
        Sniffs2Use{1} = ParsedSniffs{whichstim+1};
        [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
        plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
        set(gca,'XTick',[],'YTick',[]);
    end
    
    % plot the 5 mega odors - skip the blank
    for odor = 1:(numel(StimSettings.megaOdors)-1)
        whichsubplots = [13 19] + odor;
        subplot(nRows,nCols,whichsubplots);
        
        whichstim = find(StimulusList == StimSettings.megaOdors(odor));
        Sniffs2Use{1} = ParsedSniffs{whichstim+1};
        [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
        plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
        set(gca,'XTick',[],'YTick',[]);
    end

    set(gcf,'Position',[597   234   966   704]);
    if savefigs
        set(gcf,'Color','w');
        set(gcf,'renderer','Painters');
%         figPosition = [0.3104    0.2157    0.5031    0.7370];
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', figPosition);
        print(fullfile(myKsDir,'OdorMaps',['SniffOdorSummary',num2str(1000+n),'.eps']),'-depsc','-tiff','-r300','-painters');
%         if n == 1
%             print(fullfile(myKsDir,'OdorMaps','SniffOdorSummary.eps'),'-depsc','-tiff','-r300','-painters');
% %             exportgraphics(gcf, ...
% %                 fullfile(myKsDir,'OdorMaps','SniffOdorSummary.pdf'),...
% %                 'ContentType','vector');
%         else
%             print(fullfile(myKsDir,'OdorMaps','SniffOdorSummary.eps'),'-depsc','-tiff','-r300','-painters','-append');
% %             exportgraphics(gcf, ...
% %                 fullfile(myKsDir,'OdorMaps','SniffOdorSummary.pdf'),...
% %                 'ContentType','vector','Append',true);
%         end
        close(gcf);
    end
end

end