function [baseline,kernels,locationcoef] = ParseGaborSniffKernels(ConcatenatedKernels,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('independentcoeffs', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('kernellength', 500, @(x) isnumeric(x)); % in ms
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms

% extract values from the inputParser
params.parse(varargin{:});
notcommon       = params.Results.independentcoeffs; % independent coeffs for different odors
binsize         = params.Results.binsize;
kernellength    = floor(params.Results.kernellength/binsize); % in bins

baseline = ConcatenatedKernels(1);
% crop out the kernels
if notcommon
    locationcoef = ConcatenatedKernels(end-2:end);
else
    locationcoef = ConcatenatedKernels(end);
end

kernelparams = reshape(ConcatenatedKernels(2:end-numel(locationcoef)),[],5);
for n = 1:size(kernelparams)
    kernels(n,:) = kernelGEN_gabor(kernelparams(:,n),kernellength);
end

end