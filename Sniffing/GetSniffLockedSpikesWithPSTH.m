function [SpikeRaster, PSTHs, SniffPSTHCorr] = GetSniffLockedSpikesWithPSTH(GroupedSniffs, thisUnitSpikes, varargin)

%% extract inputs
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('yoffset', 0, @(x) isnumeric(x));
params.addParameter('window', [-0.1 1], @(x) isnumeric(x));
params.addParameter('PSTHBinsize', 20, @(x) isnumeric(x));
params.addParameter('onlyCL', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
yoffset = params.Results.yoffset;
binsize = params.Results.PSTHBinsize/1000;
window = params.Results.window;
myBins = window(1):binsize:window(2);
onlyClosedLoop = params.Results.onlyCL;


maxsniffs = 0;
for snifftype = 1:size(GroupedSniffs,2)
    SpikesPlot = [];
    SniffTS = GroupedSniffs{snifftype};
    if onlyClosedLoop
        SniffTS(SniffTS(:,8)>1,:) = [];
    end
    myPSTH = nan(size(SniffTS,1),numel(myBins));

    % plot spikes
    for x = 1:size(SniffTS,1)
        ts = SniffTS(x,1) + window;
        whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(2)));
        thisSniffSpikes = thisUnitSpikes(whichSpikes) - SniffTS(x,1);
        y = x + yoffset;
        SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes y*ones(numel(thisSniffSpikes),1)]);
        
        % for the PSTH
        thisSniffBins = -window(1) + (SniffTS(x,13) - SniffTS(x,11))*0.002; % duration in seconds
        thisSniffBins = floor(thisSniffBins/binsize);
        binMax = min(numel(myBins),thisSniffBins);
        for thisBin = 1:binMax-1
            myPSTH(x,thisBin) = numel(find(thisSniffSpikes>=myBins(thisBin) & thisSniffSpikes<myBins(thisBin+1)));
        end
    end

    %plot(SpikesPlot(:,1), SpikesPlot(:,2), '.k','Markersize', 0.5);
    SpikeRaster{snifftype} = SpikesPlot;
    maxsniffs = max(maxsniffs,x);

    % Make the average PSTHs
    % full
    PSTHs{snifftype}.mean(1,:) = mean(myPSTH,1,'omitnan') / binsize;
    nBins = sqrt(sum(~isnan(myPSTH),1));
    PSTHs{snifftype}.sem(1,:) = std(myPSTH,1,'omitnan')./nBins;
    
    if ~onlyClosedLoop
        % for each chunk
        if snifftype == 1
            maxchunks = max(unique(SniffTS(:,8)));
        end

        for chunk = 1:maxchunks
            PSTHs{snifftype}.mean(1+chunk,:)    = mean(myPSTH(find(SniffTS(:,8)==chunk),:),1,'omitnan') / binsize;
            nBins = sqrt(sum(~isnan(myPSTH(find(SniffTS(:,8)==chunk),:)),1));
            PSTHs{snifftype}.sem(1+chunk,:)     = std(myPSTH(find(SniffTS(:,8)==chunk),:),1,'omitnan')./nBins;
        end
    end

end

% get the PSTH correlations
AllStimPSTH = [];
for snifftype = 1:5
    % string together all stim types
    AllStimPSTH = horzcat(AllStimPSTH, PSTHs{snifftype}.mean(:,1:end-1)); 
end

SniffPSTHCorr = corr(AllStimPSTH');

end