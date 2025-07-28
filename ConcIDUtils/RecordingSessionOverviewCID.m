function [popFR] = RecordingSessionOverviewCID(myKsDir,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('rastermode', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('spikeheight', 0.8, @(x) isnumeric(x) && x>0 && x<=1);
params.addParameter('sessionlength', 100, @(x) isnumeric(x));
params.addParameter('FRnormalizer', 20, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
raster = params.Results.rastermode;
SpikeHeight = params.Results.spikeheight;
PSTHWindow = [0 1000*params.Results.sessionlength];
divideFRby = params.Results.FRnormalizer;
downsamplefactor = 100;

addpath(genpath('/opt/RastermapOld'));

if ~raster
    popFR = zeros(PSTHWindow(2)*downsamplefactor/1000,1);
else
    popFR = [];
end

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

%% Lets make some rasters
NumUnits = 125;
TS = 1000*[3000 3300]; % 300 seconds = 300,000 bins
timeBins = TS(1):TS(2);
myRaster = zeros(NumUnits,numel(timeBins));
TS = TS/1000;

for n = 1:NumUnits
    n = SortedBySpikes(i);

    % align to the specified event
    f = find(SingleUnits(n).spikes>TS(2),1,'first');
    if ~isempty(f)
        f = f-1;
        thisUnitSpikes = SingleUnits(n).spikes(1:f)' - TS(1);
    else
        thisUnitSpikes = SingleUnits(n).spikes' - TS(1);
    end

    %PlotRaster(thisUnitSpikes(thisUnitSpikes>0),n,[0 0 0],0.8);
    % convert spike times to milliseconds and floor values
	thisUnitSpikes = floor(1000*thisUnitSpikes);
    % remove NaNs
    thisUnitSpikes(isnan(thisUnitSpikes)) = [];
	% Make raster
    [C,~,ic] = unique(thisUnitSpikes);
    bin_counts = accumarray(ic,1);
    if ~isempty(C)
        % ignore any -ve time bins
        bin_counts((C<=0),:) = [];
        C(C<=0) = [];
        myRaster(n,C) = bin_counts;
    end
end

%% run rastermap
ops.nCall = [30 NumUnits];
ops.iPC = 1:NumUnits;
[isort1, isort2, Sm] = mapTmap(myRaster, ops);

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

for i = 1:1:size(SingleUnits,2)
   % whichunit = SortedByTetrodes(i);
    whichunit = SortedBySpikes(i);
    thistetrode = SingleUnits(whichunit).tetrode;
    if thistetrode~=whichtetrode
        MyColors = circshift(MyColors,1);
        whichtetrode = thistetrode;
    end
    if raster
        PlotRaster(SingleUnits(whichunit).spikes,i,[0 0 0],SpikeHeight);
        %PlotRaster(SingleUnits(whichunit).spikes,i,MyColors(1,:),SpikeHeight);
    else
        FR = MakePSTH(SingleUnits(whichunit).spikes', 0, PSTHWindow, 'downsample', downsamplefactor);
        veclength = min(size(popFR,1),numel(FR));
        popFR(1:veclength,1) = popFR(1:veclength,1) + FR(1:veclength,1);
        %popFR(1:veclength,2) = popFR(1:veclength,2) + FR(1:veclength,1)/max(FR);
        
        FR = i + FR/divideFRby;
        plot((1:numel(FR))/downsamplefactor,FR,'color',MyColors(1,:));
    end
end

if ~raster
    popFR = popFR/i;
    %popFR(:,1) = popFR(:,1)/divideFRby;
    %popFR(:,2) = popFR(:,2)*max(popFR(:,1));
end

end