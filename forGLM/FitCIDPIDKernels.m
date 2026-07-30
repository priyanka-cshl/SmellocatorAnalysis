function [Fitted,resnorm,residual,exitflag,output] = FitCIDPIDKernels(PIDdir, varargin)
% function to fit a kernel that when convolved with valve waveform
% recreates the PID waveform for each odor
% Inputs: Folder with the processed PID data
% from averagedPID.mat load PIDTraces (n x t), n = odors, t = samples @500hz
% and window - time before and after odor off
% from quickprocessOdorTTLs.mat get stimulus timing info and create odor
% valve waveform
    
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms

% extract values from the inputParser
params.parse(varargin{:});
binsize = params.Results.binsize;

PID     = load(fullfile(PIDdir,'quickprocessPID.mat'));
Traces  = load(fullfile(PIDdir,'averagedPID.mat'));

OdorDuration = PID.StimSettings.timing(3)/1000;
% confirm
if any(abs(diff((PID.TTLs.Trial(:,7:8))')-OdorDuration)>0.01)
    keyboard; % odor duration seems off
end

% Make Valve waveform
window = Traces.window;
timestep = binsize/1000; % in seconds
TimeVector = window(1):timestep:window(2);
ValveVector = TimeVector*0;
ValveOn = find( (TimeVector>=0) & (TimeVector<=OdorDuration) );
ValveVector(ValveOn) = 1;

% interpolate PID traces at the same resolution
sampMax = min([size(Traces.AverageOut,2) size(Traces.TimeOut,1)]);
% smoothen the PID traces to remove high freq fluctuations
PIDRaw = Traces.AverageOut(:,1:sampMax)';
fband = [0.001 200];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(500/2));
PIDfiltered = filtfilt(b,a,PIDRaw);
nOdors = size(PIDfiltered,2);
plotTraces = 1;

for i = 1:nOdors
    PIDfiltered(:,i) = PIDfiltered(:,i) - mean(PIDfiltered(1:500*abs(window(1)),i));
    PIDfiltered(:,i) = PIDfiltered(:,i)/max(PIDfiltered(:,i));
end

[PIDtraces] = interp1(Traces.TimeOut(1:sampMax,1),PIDfiltered,TimeVector);
PIDtraces(1,:) = [];
TimeVector(:,1) = [];
ValveVector(:,1) = [];

if plotTraces
   figure; 
    for i = 1:16
        subplot(4,4,i);
        plot(Traces.TimeOut(1:sampMax,1),PIDfiltered(:,i),'b');
        hold on
        plot(TimeVector,PIDtraces(:,i),'k');
        set(gca,'XLim',[-0.5 2.5]);
    end
end
nOdors = size(PIDtraces,2);
xdata = repmat(ValveVector',1,nOdors);
ydata = PIDtraces;

% start clock
tic

%% 1 : set up the error minimization
StartingKernels = zeros(5/timestep,nOdors);
model_fit = @pid_out; 
lb = -100 + 0*StartingKernels;
ub = 100 + 0*StartingKernels;
Eval_max = 1e+6; Iter_max = 1e+6;  
Fun_Tol  = 1e-8; Step_Tol = 1e-8;


options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max); %,'TolFun',Fun_Tol,'TolX',Step_Tol);
%options = optimset('Algorithm','levenberg-marquardt');
%[kernelsout] = lsqcurvefit(model_fit,StartingKernels,xdata,ydata,lb,ub,options);
[kernelsout,resnorm,residual,exitflag,output] = lsqcurvefit(model_fit,StartingKernels,xdata,ydata,lb,ub,options);

%% 2 : define the convolution function (model_fit)
    function [zdata] = pid_out(StartingKernels,xdata)

        % make loal copies
        Starting_Kernels    = StartingKernels;
        x_data              = xdata;
        maxT                = size(x_data,1)+1;
        convmode = 'full'; %'same';
        for k = 1:size(x_data,2) 
            zdata(:,k) = conv(x_data(:,k)',Starting_Kernels(:,k)',convmode);
        end

        zdata(maxT:end,:) = [];
    end

%% stop clock
toc

end





