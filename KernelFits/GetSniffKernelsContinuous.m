function [kernelsout] = GetSniffKernelsContinuous(StartingKernels, Predictors, ObservedPSTHs, varargin)
% function to the best kernel that explains the data

% StartingKernels
    % = [baseline ITIKernel AirKernel OdorKernel locationcoef];
    
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms
params.addParameter('rectifyFR', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
rectifyFR = params.Results.rectifyFR;

% start clock
tic

xdata = Predictors;
ydata = ObservedPSTHs;

%% 1 : set up the error minimization
Eval_max = 1e+6; Iter_max = 1e+6; model_fit = @sniff_out; lb = -100 + 0*StartingKernels; ub = 100 + 0*StartingKernels;
Fun_Tol  = 1e-8; Step_Tol = 1e-8;
options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max); %,'TolFun',Fun_Tol,'TolX',Step_Tol);
%options = optimset('Algorithm','levenberg-marquardt');
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
    [baseline,kernels,locationcoef] = ParseSniffKernels(Starting_Kernels);
    % get psth
    [zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,x_data);
    
    if rectifyFR
        zdata(zdata<0) = 0;
    end
end
%zdata = [];

%% stop clock
toc

end





