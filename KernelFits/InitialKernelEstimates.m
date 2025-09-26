function [StartingKernels] = InitialKernelEstimates(SniffPSTHs, SniffParams, varargin)
% function to get starting kernel conditions from the sniffaligned PSTHs

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('kernellength', 500, @(x) isnumeric(x)); % in ms
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms
params.addParameter('independentcoeffs', 0, @(x) logical(x));

% extract values from the inputParser
params.parse(varargin{:});
binsize         = params.Results.binsize;
kernellength    = floor(params.Results.kernellength/binsize); % in bins
notcommon = params.Results.independentcoeffs; % independent coeffs for different odors

%% Step 1: extract the sniff length column 
%   and convert padded columns to NaNs

Snifflengths = SniffPSTHs(:,1) + 1;
SniffPSTHs(:,1) = [];
for i = 1:size(SniffPSTHs,1) % every psth
    SniffPSTHs(i,Snifflengths(i):end) = NaN;
end

%% Step 2: estimate ITI, air, odor kernels and baseline
ITISniffs       = find(SniffParams(:,1)==-1);
ITIKernel       = mean(SniffPSTHs(ITISniffs,:),'omitnan');
ITIKernel(isnan(ITIKernel)) = 0;
baseline        = mean(ITIKernel((kernellength+1):end),'omitnan');
if isnan(baseline)
    baseline = 0;
end

FarLocations    = 60; % floor(prctile(SniffParams(find(SniffParams(:,1)>0),2),95))
AirSniffs       = intersect(find(SniffParams(:,1)>=0), find(abs(SniffParams(:,2))>FarLocations));
AirKernel       = mean(SniffPSTHs(AirSniffs,:),'omitnan');
AirKernel(isnan(AirKernel)) = 0;

NearLocations   = 8;
for odor = 1:3
    OdorSniffs      = intersect(find(SniffParams(:,1)==odor), find(abs(SniffParams(:,2))<=NearLocations));
    OdorKernel{odor}      = mean(SniffPSTHs(OdorSniffs,:),'omitnan');
    OdorKernel{odor}(isnan(OdorKernel{odor})) = 0;
end
locationcoef    = -0.08; % C(x) = exp(-0.1729/2*location(a.u.))

% crop kernels 
ITIKernel(:,(kernellength+1):end) = [];
ITIKernel = ITIKernel - baseline;
AirKernel(:,(kernellength+1):end) = [];
if any(AirKernel)
    AirKernel = AirKernel - ITIKernel - baseline; % effect of Air alone?
end
for odor = 1:3
    OdorKernel{odor}(:,(kernellength+1):end) = [];
    if any(OdorKernel{odor})
        OdorKernel{odor} = OdorKernel{odor} - ITIKernel - AirKernel - baseline; % effect of odor alone
    end
end
if notcommon
    StartingKernels = [baseline ITIKernel AirKernel cell2mat(OdorKernel) locationcoef locationcoef locationcoef]; % independent coeffs for each odor
else
    StartingKernels = [baseline ITIKernel AirKernel cell2mat(OdorKernel) locationcoef];
end






