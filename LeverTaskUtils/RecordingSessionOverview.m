function [popFR] = RecordingSessionOverview(SingleUnits,varargin)

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

if ~raster
    popFR = zeros(PSTHWindow(2)*downsamplefactor/1000,2);
else
    popFR = [];
end

for i = 1:size(SingleUnits,2)
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id];
end

[~,SortedByTetrodes] = sort(foo(:,1));
ClusterIds = foo(SortedByTetrodes,2);
whichtetrode = SingleUnits(SortedByTetrodes(1)).tetrode;
MyColors(1,:) = Plot_Colors('r');
MyColors(2,:) = Plot_Colors('k');

for i = 1:1:size(SingleUnits,2)
    whichunit = SortedByTetrodes(i);
    thistetrode = SingleUnits(whichunit).tetrode;
    if thistetrode~=whichtetrode
        MyColors = circshift(MyColors,1);
        whichtetrode = thistetrode;
    end
    if raster
        PlotRaster(SingleUnits(whichunit).spikes,i,MyColors(1,:),SpikeHeight);
    else
        FR = MakePSTH(SingleUnits(whichunit).spikes', 0, PSTHWindow, 'downsample', downsamplefactor);
        veclength = min(size(popFR,1),numel(FR));
        popFR(1:veclength,1) = popFR(1:veclength,1) + FR(1:veclength,1);
        popFR(1:veclength,2) = popFR(1:veclength,2) + FR(1:veclength,1)/max(FR);
        
        FR = i + FR/divideFRby;
        plot((1:numel(FR))/downsamplefactor,FR,'color',MyColors(1,:));
    end
end

if ~raster
    popFR = popFR/i;
    %popFR(:,1) = popFR(:,1)/divideFRby;
    popFR(:,2) = popFR(:,2)*max(popFR(:,1));
end
end