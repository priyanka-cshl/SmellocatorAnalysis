%% script to extract the sniff aligned spiking responses
%  of a given neuron 

%% Step 1: Data Paths
myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

%% Get the sniff timestamps
%  Load Sniffs from the T3 session and make a vector with binarized sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
load(fullfile(myKsDir,'SessionDetails.mat'));

if exist('CuratedSniffTimestamps','var')
    AllSniffs = CuratedSniffTimestamps(find(CuratedSniffTimestamps(:,8)>0),:);
    AllSniffs(:,4:7) = 0;
    AllSniffs(:,3) = AllSniffs(:,3) - AllSniffs(:,1);
end

% add odor TTL info if available
if exist(fullfile(myKsDir,'quickprocessOdorTTLs.mat'))
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
    
end

startTS = 0;
for i = 1:numel(Files.Samples)
    endTS = Files.Samples(i)/30000 + startTS; % in seconds
    whichsniffs = find( (AllSniffs(:,1)>=startTS) & (AllSniffs(:,1)<endTS) );
    AllSniffs(whichsniffs,8) = i; % session phase
    startTS = endTS;
end

% group sniffs by odors
% there are no air OFF sniffs here
AllSniffs(:,4) = 1; % air is always on
[ParsedSniffs, StimulusList] = ParseSniffsByStimuli(AllSniffs, 'SortBy', 3);
Sniffs2Use{1} = ParsedSniffs{2}; % Air sniffs
% % Sniffs2Use{1}(:,8) = 1;
% Sniffs2Use{1}(11000:end,:) = [];
% Sniffs2Use{1}(5000:8500,:) = [];

%% Get the PSTHs at the same resolution
%  Load spikes
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

%selectedUnits = [19 24 29 30 54 56 57 101 108 133]; % 134 167
selectedUnits = [19 54 24 56 57 29 30 108 101 133];
selectedUnits = [30 133 54 56 101];


%%
savefigs = 1;
nRows = 6*3; nCols = 10; 
nTOT = numel(selectedUnits); % nUnits %

for thisunit = 1:nTOT
    
    if ~rem(thisunit-1,nCols);
        figure;
    end

    % actual spikes
    if nTOT == nUnits
        k = thisunit;
    else
        k = selectedUnits(thisunit);
    end
    
    thisUnitSpikes = SingleUnits(k).spikes;
    % sort by duration only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},3,'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    whichplots = (1:10:(nRows*nCols/3)) + rem(thisunit-1,nCols);
    subplot(nRows+1,nCols,whichplots);
    %plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    h1 = plot(nan,nan,'.k','Markersize', 0.5);
    h1.XData = SpikeRaster{1}(:,1)'; 
    h1.YData = SpikeRaster{1}(:,2)';
    hold on
    plot(Sniffs2Use{1}(:,1)*0,1:size(Sniffs2Use{1},1),'r');
    plot(Sniffs2Use{1}(:,3),1:size(Sniffs2Use{1},1),'r');

    % sort by session phase and duration only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},[8 3],'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);

    whichplots = whichplots + (nRows*nCols/3);
    subplot(nRows+1,nCols,whichplots);
    %plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    h2 = plot(nan,nan,'.k','Markersize', 0.5);
    h2.XData = SpikeRaster{1}(:,1)';
    h2.YData = SpikeRaster{1}(:,2)';
    hold on
    plot(Sniffs2Use{1}(:,1)*0,1:size(Sniffs2Use{1},1),'r');
    plot(Sniffs2Use{1}(:,3),1:size(Sniffs2Use{1},1),'r');

    % sort by session phase and occurence only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},[8 1],'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);

    whichplots = whichplots + (nRows*nCols/3);
    subplot(nRows+1,nCols,whichplots);
    %plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    h3 = plot(nan,nan,'.k','Markersize', 0.5);
    h3.XData = SpikeRaster{1}(:,1)';
    h3.YData = SpikeRaster{1}(:,2)';
%     hold on
%     plot(Sniffs2Use{1}(:,1)*0,1:size(Sniffs2Use{1},1),'r');
%     plot(Sniffs2Use{1}(:,3),1:size(Sniffs2Use{1},1),'r');

    if ~rem(thisunit,nCols) || thisunit == nTOT
        set(gcf,'Position',[2015 6 1748 989]);
        if savefigs
            set(gcf,'Color','w');
            set(gcf,'renderer','Painters');
            print(['/mnt/data/simulation/NewSelectedObservedAirResponse',num2str(thisunit),'.eps'],'-depsc','-tiff','-r300','-painters');
            close(gcf);
        end
    end
    
    
end
