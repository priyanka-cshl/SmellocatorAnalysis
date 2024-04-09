function [kernelsout,resnorm,residual,exitflag,output] = GetSniffKernels(StartingKernels, SniffParams, SniffPSTHs, varargin)
% function to the best kernel that explains the data

% StartingKernels
    % = [baseline ITIKernel AirKernel OdorKernel locationcoef];
    
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms


% extract values from the inputParser
params.parse(varargin{:});
binsize = params.Results.binsize;

% start clock
tic

xdata = [SniffParams SniffPSTHs(:,1) floor(abs(SniffParams(:,9))*1000/binsize)];
ydata = SniffPSTHs(:,2:end);

%% 1 : set up the error minimization
Eval_max = 1e+6; Iter_max = 1e+6; model_fit = @sniff_out; lb = []; ub = [];
%options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max);
options = optimset('Algorithm','levenberg-marquardt');
%[kernelsout] = lsqcurvefit(model_fit,StartingKernels,xdata,ydata,lb,ub,options);
[kernelsout,resnorm,residual,exitflag,output] = lsqcurvefit(model_fit,StartingKernels,xdata,ydata,lb,ub,options);

% sniff params
    % 1-4: [currsniffstate currsniffloc currinhend currsniffend ...
    % 5-6:  currsniffTrialID currsniffIndex ...
    % 7:10: prevsniffstate prevsniffloc previnhstart previnhend ...
    % 11:   Snifflength in bins ...
    % 12:   PrevSniffLength in bins]
    
%% 2 : define the convolution function (model_fit)
function [zdata] = sniff_out(StartingKernels,xdata)
    
    % make loal copies
    Starting_Kernels    = StartingKernels;
    x_data              = xdata;
    
    % crop out the kernels 
    baseline        = Starting_Kernels(1);
    locationcoef    = Starting_Kernels(end);
    kernels         = reshape(Starting_Kernels(2:end-1),[],5);
    
    snifflengths    = x_data(:,11:12);
    
    z_data = zeros(size(x_data,1),max(x_data(:,11)));
    for i = 1:size(x_data,1) % every sniff
        % this sniff
        stimstate1 = x_data(i,1); % -1 for ITI, 0, 1, 2, 3 for air, odor 1-3
        location1  = exp(abs(x_data(i,2)) * locationcoef);
        
        R1 = kernels(:,1) + ...
            (stimstate1>=0)*kernels(:,2) + ...
            (stimstate1==1)*location1*kernels(:,3) + ...
            (stimstate1==2)*location1*kernels(:,4) + ...
            (stimstate1==3)*location1*kernels(:,5);
        
        % previous sniff
        stimstate2 = x_data(i,7); % -1 for ITI, 0, 1, 2, 3 for air, odor 1-3
        location2  = exp(abs(x_data(i,8)) * locationcoef);
        
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
        z_data(i,1:length(R1)) = R1;        
    end
    
    zdata = z_data;
    zdata(zdata<0) = 0;
    for k = 1:size(zdata,1)
        zdata(k,(1+snifflengths(k,1)):end) = 0;
    end
end
zdata = [];

%% stop clock
toc

end





