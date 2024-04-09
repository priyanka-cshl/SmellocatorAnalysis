function [PSTHOut] = GetPredictedSniffPSTH(SniffParams,EstimatedKernels,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms


% extract values from the inputParser
params.parse(varargin{:});
binsize = params.Results.binsize;

% start clock
tic

SniffParams = [SniffParams floor(abs(SniffParams(:,4))*1000/binsize) floor(abs(SniffParams(:,9))*1000/binsize)];

% crop out the kernels
baseline        = EstimatedKernels(1);
locationcoef    = EstimatedKernels(end);
kernels         = reshape(EstimatedKernels(2:end-1),[],5);

snifflengths    = SniffParams(:,11:12);

PSTHOut = zeros(size(SniffParams,1),max(SniffParams(:,11)));
for i = 1:size(SniffParams,1) % every sniff
    % this sniff
    stimstate1 = SniffParams(i,1); % -1 for ITI, 0, 1, 2, 3 for air, odor 1-3
    location1  = exp(abs(SniffParams(i,2)) * locationcoef);
    
    R1 = kernels(:,1) + ...
        (stimstate1>=0)*kernels(:,2) + ...
        (stimstate1==1)*location1*kernels(:,3) + ...
        (stimstate1==2)*location1*kernels(:,4) + ...
        (stimstate1==3)*location1*kernels(:,5);
    
    % previous sniff
    stimstate2 = SniffParams(i,7); % -1 for ITI, 0, 1, 2, 3 for air, odor 1-3
    location2  = exp(abs(SniffParams(i,8)) * locationcoef);
    
    R2 = kernels(:,1) + ...
        (stimstate2>=0)*kernels(:,2) + ...
        (stimstate2==1)*location2*kernels(:,3) + ...
        (stimstate2==2)*location2*kernels(:,4) + ...
        (stimstate2==3)*location2*kernels(:,5);
    
    bins2delete = snifflengths(i,2);
    if length(R2)>bins2delete
        R2(1:bins2delete,:) = [];
        R2(end+1:length(R1),:) = 0;
    else
        R2 = R2*0;
    end
    
    R1 = R1 + R2 + baseline;
    PSTHOut(i,1:length(R1)) = R1;
end

PSTHOut(PSTHOut<0) = 0;
for k = 1:size(PSTHOut,1)
    PSTHOut(k,(1+snifflengths(k,1)):end) = 0;
end
end