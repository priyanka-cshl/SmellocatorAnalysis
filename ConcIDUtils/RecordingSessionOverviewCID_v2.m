% function [popFR] = RecordingSessionOverviewCID(myKsDir,varargin)
% 
% %% parse input arguments
% narginchk(1,inf)
% params = inputParser;
% params.CaseSensitive = false;
% params.addParameter('rastermode', true, @(x) islogical(x) || x==0 || x==1);
% params.addParameter('spikeheight', 0.8, @(x) isnumeric(x) && x>0 && x<=1);
% params.addParameter('sessionlength', 100, @(x) isnumeric(x));
% params.addParameter('FRnormalizer', 20, @(x) isnumeric(x));
% 
% % extract values from the inputParser
% params.parse(varargin{:});
% raster = params.Results.rastermode;
% SpikeHeight = params.Results.spikeheight;
% PSTHWindow = [0 1000*params.Results.sessionlength];
% divideFRby = params.Results.FRnormalizer;
% downsamplefactor = 100;

myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

raster = 1; 
SpikeHeight = 0.8;
PSTHWindow = [0 1000*100];
divideFRby = 20;
downsamplefactor = 100;


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

%% useful?
for i = 1:size(SingleUnits,2)
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id numel(SingleUnits(i).spikes)];
end

[~,SortedByTetrodes] = sort(foo(:,1));
ClusterIds = foo(SortedByTetrodes,2);
whichtetrode = SingleUnits(SortedByTetrodes(1)).tetrode;
MyColors(1,:) = Plot_Colors('r');
MyColors(2,:) = Plot_Colors('k');

[~,SortedBySpikes] = sort(foo(:,3));
ClusterIds = foo(SortedBySpikes,2);

%% actual plotting?
BoxColors = brewermap(20, 'Set3');

%%
% plot odor boxes?
OdorList = unique(AllSniffs(:,5));
OdorList(find(~OdorList),:) = [];
figure;
hold on;
for odorct = 1:numel(OdorList)
    odor = OdorList(odorct);
    %h.(['OdorPlot',num2str(odor)]) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(odorct)]));
    h.(['OdorPlot',num2str(odor)]) = fill(NaN,NaN,BoxColors(odorct,:));
    h.(['OdorPlot',num2str(odor)]).EdgeColor = 'none';
    thisOdorTrials = TTLs.Trial(find(TTLs.Trial(:,4)==odor),:);
    ValveTS = thisOdorTrials(:,7:8)';
    h.(['OdorPlot',num2str(odor)]).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat((10+nUnits)*[0 1 1 0]',size(ValveTS,2),1)];
    h.(['OdorPlot',num2str(odor)]).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';

end

% plot sniff boxes
h.sniffplot = fill(NaN,NaN,BoxColors(odorct+1,:));
h.sniffplot.EdgeColor = 'none';
h.sniffplot.FaceAlpha = 0.5;
ValveTS = AllSniffs(:,1:2)';
h.sniffplot.Vertices = [ ...
    reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
    repmat((10+nUnits)*[0 1 1 0]',size(ValveTS,2),1)];
h.sniffplot.Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';

%end
% order according to specific sniff
% last sniff before 542
whichsniff = find(AllSniffs(:,1)<542,1,'last');
whichsniff = find(AllSniffs(:,1)<868,1,'last');

whichsniff = find(AllSniffs(:,1)<1124,2,'last');
% whichsniff = find(AllSniffs(:,1)>706,1,'first');
whichsniff(1:5,1) = whichsniff(1):(whichsniff(1)+4);

% whichsniff = find(AllSniffs(:,1)>740,2,'first');

for i = 1:1:size(SingleUnits,2)
    for m = 1:size(whichsniff,1)
        sniffstart = AllSniffs(whichsniff(m),1);
        firstspike =  find(SingleUnits(i).spikes>=sniffstart,1,'first');
        if ~isempty(firstspike)
            delta(i,m) = SingleUnits(i).spikes(firstspike) - sniffstart;
        else
            delta(i,m) = -1;
        end
    end
end
[~,SortedByLatency] = sort(mean(delta,2),'descend');

for i = 1:1:size(SingleUnits,2)
    whichunit = SortedByLatency(i);
    PlotRaster(SingleUnits(whichunit).spikes,i,[0 0 0],SpikeHeight);
end

%%
%%
for i = 1:1:size(SingleUnits,2)
    % whichunit = SortedByTetrodes(i);
    whichunit = SortedBySpikes(i);
    thistetrode = SingleUnits(whichunit).tetrode;
    PlotRaster(SingleUnits(whichunit).spikes,i,[0 0 0],SpikeHeight);
    %PlotRaster(SingleUnits(whichunit).spikes,i,MyColors(1,:),SpikeHeight);

end