function [] = RecordingSessionOverview(SingleUnits,varargin)

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
            FR = MakePSTH(SingleUnits(whichunit).spikes', 0, PSTHWindow);
            FR = i + FR/divideFRby;
            plot((1:numel(FR))/1000,FR,'color',MyColors(1,:));
        end
    end
end