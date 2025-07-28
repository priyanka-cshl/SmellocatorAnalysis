
myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

raster = 1; 
SpikeHeight = 1; %0.8;
BoxColors = brewermap(20, 'Set3');

%% load stuff
% load single units
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

% load Sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
load(fullfile(myKsDir,'SessionDetails.mat'));
if exist('CuratedSniffTimestamps','var')
    AllSniffs = CuratedSniffTimestamps(find(CuratedSniffTimestamps(:,8)>0),:);
    AllSniffs(:,4:7) = 0;
    AllSniffs(:,3) = AllSniffs(:,3) - AllSniffs(:,1);
end

% load odor TTL info if available
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

OdorList = unique(AllSniffs(:,5));
OdorList(find(~OdorList),:) = [];

%% useful?
for i = 1:nUnits
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id numel(SingleUnits(i).spikes)];
end
[~,SortedByTetrodes] = sort(foo(:,1));
[~,SortedBySpikes] = sort(foo(:,3));

%% for extracting time windows around odor periods
myOdor = 15;
myOdorTrials = find(TTLs.Trial(:,4)==myOdor);
myOdorStarts = (TTLs.Trial(myOdorTrials,7));
myOdorStops = (TTLs.Trial(myOdorTrials,8));

% Sort by aligning to a specific sniff
% % last sniff before a given TS = 542
% TS = 542; % 868
% whichsniff = find(AllSniffs(:,1)<TS,1,'last');

% % last n sniffs before a given TS = 1124
% TS = 1124; n = 2; nT = 5; % take the sequence of m sniffs after the desired sniff
% %TS = 612; n = 4; nT = 5; % take the sequence of m sniffs after the desired sniff
% whichsniff = find(AllSniffs(:,1)<TS,n,'last');
% whichsniff(1:(nT+1),1) = whichsniff(1):(whichsniff(1)+nT); 

TS = myOdorStops(2); n = 6; nT = 5; % take the sequence of m sniffs after the desired sniff
whichsniff = find(AllSniffs(:,1)<TS,n,'last');

for i = 1:nUnits
    for m = 1:nT
        sniffstart = AllSniffs(whichsniff(m),1); 
        sniffend   = AllSniffs(whichsniff(m+1),1); 
        firstspike = find(SingleUnits(i).spikes>=sniffstart,1,'first');
        lastspike  = find(SingleUnits(i).spikes>=sniffend,1,'first');
        if ~isempty(firstspike)
            spikecount(i,m) = lastspike - firstspike;
            delta(i,m) = SingleUnits(i).spikes(firstspike) - sniffstart;
        else
            spikecount(i,m) = -1;
            delta(i,m) = -1;
        end
    end
end
M = [mean(delta,2) mean(spikecount,2)];
[~,SortOrder] = sortrows(M,[2 1],'ascend');
%[~,SortOrder] = sortrows(M,[1 2],'descend');
%[~,SortedByLatency] = sort(mean(delta,2),'descend');

% actual plotting
tTrials = 3;
whichTrials = [1 2 3];
figure;
yLims = [88 153]; 
yLims = [68 158];

for t = 1:tTrials
    subplot(tTrials,1,t); hold on;

    % plot odor boxes
    for odorct = 1:numel(OdorList)
        odor = OdorList(odorct);
        thisOdorTrials = TTLs.Trial(find(TTLs.Trial(:,4)==odor),:);
        ValveTS = thisOdorTrials(:,7:8)';
        plotEventBoxes(['OdorPlot',num2str(odor)],ValveTS,BoxColors(odorct,:),1,nUnits);
    end

    % plot sniff boxes
    plotEventBoxes('sniffplot',AllSniffs(:,1:2)',BoxColors(odorct+1,:),0.5,nUnits);
    
    % plot spikes
    for i =  73:nUnits %88:nUnits
        %whichunit = SortedBySpikes(i);
        whichunit = SortOrder(i);
        %PlotRaster(SingleUnits(whichunit).spikes,i,[0 0 0],SpikeHeight);
        PlotRaster_v2(SingleUnits(whichunit).spikes,i,[0 0 0],'ticklength',SpikeHeight,'tickwidth',0.5);
    end

    set(gca,'XLim', myOdorStarts(whichTrials(t)) + [-4 8],...
        'YLim', yLims, 'TickDir', 'out');
    set(gcf,'Position',[1990 50 510 900]);
end


%% functions
function [] = plotEventBoxes(plotHandle,ValveTS,BoxColor,Transparency,BoxHeight)
    
    h.(plotHandle) = fill(NaN,NaN,BoxColor);
    h.(plotHandle).EdgeColor = 'none';
    h.(plotHandle).FaceAlpha = Transparency;
    h.(plotHandle).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat((10+BoxHeight)*[0 1 1 0]',size(ValveTS,2),1)];
    h.(plotHandle).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end