function [AllPSTHs,SniffParams] = SniffAlignedFRs(AllSniffs, thisUnitSpikes, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('kernelsize', 100, @(x) isnumeric(x));
params.addParameter('bufferwindow', 500, @(x) isnumeric(x));
params.addParameter('smooth', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('binsize', 50, @(x) isnumeric(x));
params.addParameter('psthtype', 'binned', @(x) any(validatestring(x, {'raw', 'binned'})));

% extract values from the inputParser
params.parse(varargin{:});
kernelsize      = params.Results.kernelsize;
bufferwindow    = params.Results.bufferwindow/1000;
smoothpsth      = params.Results.smooth;
binsize         = params.Results.binsize;
psthtype        = params.Results.psthtype;
if ~strcmp(psthtype,'binned')
    binsize = 1;
end

% for psth smoothing
taxis = (-500:500)/binsize;  % make a time axis of 1000 ms
gauss_kernel = normpdf(taxis, 0, kernelsize/binsize);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

AllPSTHs = []; SniffParams = [];

% Get Spikes for each sniffs, convert to PSTH, also get sniff params for
% this and the previous sniff
for x = 1:size(AllSniffs,1)
    ts = AllSniffs(x,[5 6 7 8 9]); % [previnhstart previnhend thisinhstart thisinhend thissniffend]
    t1 = ts(3) - bufferwindow; % for psth smoothing
    t2 = ts(5) + bufferwindow;
    whichspikes = intersect(find(thisUnitSpikes>=t1), find(thisUnitSpikes<=t2));
    whichspikes = thisUnitSpikes(whichspikes) - t1; % ref w.r.t. sniffstart and no negative spiketimes 
    snifflength = ts(5) - ts(3); % in seconds
    
    rasterlength = floor((t2-t1)*1000/binsize);
    myRaster    = zeros(1,rasterlength);
    % convert spiketimes to raster
    [C,~,ic] = unique(floor(1000*whichspikes/binsize));
    bin_counts = accumarray(ic,1);
    if ~isempty(C)
        % ignore any -ve time bins
        bin_counts((C<=0),:) = [];
        C(C<=0) = [];
        myRaster(1,C) = bin_counts; % raster at the given binsize resolution
    end
        
    if smoothpsth
        % make the psth at the desired resolution
        myFR = (1000/binsize)*conv(myRaster,gauss_kernel,'same'); % in Hz
    else
        myFR = myRaster;
    end
    
    % remove the bufferperiods
    extrabins = (bufferwindow*1000)/binsize;
    myFR(:,1:extrabins) = [];
    myFR(:,(1+end-extrabins):end) = [];
    % append snifflength in the beginning
    myFR = [floor(snifflength*1000/binsize) myFR];
    
    AllPSTHs(x,1:numel(myFR)) = myFR;
    SniffParams(x,:) = [ AllSniffs(x,[2 13]) ts(4:5) - ts(3) ...
        AllSniffs(x,[1 3]) ...
        AllSniffs(x,[18 12]) ts(1:2) - ts(3)]; 
        % currsniffstate currsniffloc currinhend currsniffend 
        % currsniffTrialID currsniffIndex
        % prevsniffstate prevsniffloc previnhstart previnhend 
end
end