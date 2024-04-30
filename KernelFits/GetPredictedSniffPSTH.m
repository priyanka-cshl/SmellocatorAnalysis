function [PSTHOut,kernels,locationcoef] = GetPredictedSniffPSTH(SniffParams,EstimatedKernels,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms

% extract values from the inputParser
params.parse(varargin{:});
binsize = params.Results.binsize;

% start clock

SniffParams = [SniffParams floor(abs(SniffParams(:,4))*1000/binsize) floor(abs(SniffParams(:,9))*1000/binsize)];
if rem(numel(EstimatedKernels)-2,5) == 0
    notcommon = 0;
else
    notcommon = 1;
end

% crop kernels
[baseline,kernels,locationcoef] = ParseSniffKernels(EstimatedKernels,'independentcoeffs',notcommon);
% get psth
PSTHOut = SniffKernels2PSTH(baseline,kernels,locationcoef,SniffParams);

end