function [SpikeRaster, maxsniffs, psth] = GetSniffLockedSpikesAndPSTH(GroupedSniffs, thisUnitSpikes, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('yoffset', 0, @(x) isnumeric(x));
params.addParameter('window', [-0.1 1], @(x) isnumeric(x));
params.addParameter('binsize', 0.01, @(x) isnumeric(x));


% extract values from the inputParser
params.parse(varargin{:});
yoffset = params.Results.yoffset;
binsize = params.Results.binsize;
window  = params.Results.window;

maxsniffs = 0;
for snifftype = 1:size(GroupedSniffs,2)
    SpikesPlot = [];
    SniffTS = GroupedSniffs{snifftype};

    % plot spikes
    for x = 1:size(SniffTS,1)
        ts = SniffTS(x,1) + window;
        whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(2)));
        thisSniffSpikes = thisUnitSpikes(whichSpikes) - SniffTS(x,1);
        y = x + yoffset;
%         if SniffTS(x,8) == 1
            SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes y*ones(numel(thisSniffSpikes),1)]);
%         else
%             SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes -(y+100*(SniffTS(x,8)-1))*ones(numel(thisSniffSpikes),1)]);
%         end
    end

    %plot(SpikesPlot(:,1), SpikesPlot(:,2), '.k','Markersize', 0.5);
    SpikeRaster{snifftype} = SpikesPlot;
    maxsniffs = max(maxsniffs,x);

    % make psth
    myPSTH = [];
    for bins = window(1):binsize:(window(2)-binsize)
        nSpikes = find((SpikesPlot(:,1)>=bins)&(SpikesPlot(:,1)<(bins+binsize)));
        myPSTH = [myPSTH, numel(nSpikes)];
    end
    myPSTH = myPSTH/maxsniffs;
    psth{snifftype} = myPSTH;
end

end