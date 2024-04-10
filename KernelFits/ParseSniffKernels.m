function [baseline,kernels,locationcoef] = ParseSniffKernels(ConcatenatedKernels,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('independentcoeffs', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
notcommon = params.Results.independentcoeffs; % independent coeffs for different odors

baseline = ConcatenatedKernels(1);
% crop out the kernels
if notcommon
    locationcoef = ConcatenatedKernels(end-2:end);
else
    locationcoef = ConcatenatedKernels(end);
end
kernels = reshape(ConcatenatedKernels(2:end-numel(locationcoef)),[],5);

end