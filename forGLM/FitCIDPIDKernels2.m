function [Fitted,resnorm,residual,exitflag,output] = FitCIDPIDKernels2(PIDdir, varargin)
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
startat = ValveOn(1);
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
PIDtraces(1:startat,:) = [];
TimeVector(:,1:startat) = [];
ValveVector(:,1:startat) = [];

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
KernelMode = 'biphasic';

switch KernelMode
    case 'full'
        StartingKernels = zeros(5/timestep,nOdors);
        model_fit = @pid_out;
        lb = -100 + 0*StartingKernels;
        ub = 100 + 0*StartingKernels;
        Eval_max = 1e+6; Iter_max = 1e+6;
        Fun_Tol  = 1e-8; Step_Tol = 1e-8;
        options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max); %,'TolFun',Fun_Tol,'TolX',Step_Tol);
        [kernelsout,resnorm,residual,exitflag,output] = lsqcurvefit(model_fit,StartingKernels,xdata,ydata,lb,ub,options);
    case 'exponentials'
        nParams = 4;
        StartingParams = repmat([1, 0.05, 0.3, 0], nOdors, 1)'; % nParams x nOdors, tune seeds
        lb = repmat([0,    0.001, 0.001, 0]',   1, nOdors);
        ub = repmat([5,    2,     5,     0.5]', 1, nOdors);
        model_fit = @pid_param_out;
        options = optimset('MaxFunEvals',1e6,'MaxIter',1e6);
        [paramsout,resnorm,residual,exitflag,output] = ...
            lsqcurvefit(model_fit, StartingParams, xdata, ydata, lb, ub, options);
        [kernelsout,PIDOut] = makeKenel(ValveVector,paramsout,timestep);
    case 'adaptive'
        % params per odor: [A0, Ass, tau_rise, tau_adapt, delay]
        nParams = 5;
        StartingParams = repmat([1, 0.7, 0.05, 1.0, 0]', 1, nOdors); % seed guesses — tune these
        lb = repmat([0,   0,   0.001, 0.01, 0]',  1, nOdors);
        ub = repmat([5,   5,   2,     10,   0.5]', 1, nOdors);
        model_fit = @pid_adapt_out;
        options = optimset('MaxFunEvals',1e6,'MaxIter',1e6);
        [paramsout,resnorm,residual,exitflag,output] = ...
            lsqcurvefit(model_fit, StartingParams, xdata, ydata, lb, ub, options);
        [kernelsout,PIDOut] = makeKenel(ValveVector,paramsout,timestep);
    case 'biphasic'
        % params per odor: [A1, tau1, n1, A2, dtau, n2, delay]
        % note: tau2 = tau1 + dtau internally (guarantees tau2 > tau1,
        % avoiding the lobe sign-flip degeneracy) WITHOUT imposing an
        % artificial absolute floor on tau2.
        nParams = 7;

        lb1 = [0,    0.01, 1, 0,    0.02, 1, 0]';
        ub1 = [5,    1,    6, 5,    6,    6, 0.5]';
        lb  = repmat(lb1, 1, nOdors);
        ub  = repmat(ub1, 1, nOdors);

        % lambda: regularization weight on A2 (the negative/adapting lobe).
        % A2 is only allowed to grow away from 0 if doing so meaningfully
        % reduces the fit residual -- this keeps odors that are genuinely
        % monophasic (no sag/undershoot) from picking up a spurious small
        % negative lobe that shows up as a dip after valve-off.
        lambda = 0.05; % tune: larger = more aggressive suppression of A2

        % Multi-start seeds: [A1, tau1, n1, A2, dtau, n2, delay]
        % Cover a range of onset/adaptation speeds so odors with fast
        % adaptation (sag starting within ~0.5s of onset), slow/broad
        % adaptation (gradual decline over ~1-1.5s), and no adaptation
        % all have a seed that's in the right basin.
        seedList = [ ...
            1,   0.10, 2, 0.30, 0.30, 2, 0;   % fast onset, fast adaptation
            1,   0.10, 2, 0.30, 1.50, 2, 0;   % fast onset, slow adaptation
            1,   0.30, 2, 0.30, 0.30, 2, 0;   % slower onset, fast adaptation
            1,   0.05, 3, 0.05, 0.50, 2, 0;   % near-monophasic (small A2)
            1,   0.15, 2, 0.70, 1.20, 3, 0;   % strong, slow+broad adaptation
            1,   0.15, 2, 1.00, 0.80, 2, 0;   % strong, moderate-speed adaptation
            1,   0.10, 3, 0.90, 2.00, 4, 0;   % strong, very slow+broad adaptation
            1,   0.10, 2, 0.25, 3.00, 3, 0;   % weak, very slow suppression (still rising within window)
            1,   0.10, 2, 0.15, 4.50, 3, 0;   % weak, near-linear suppression (tau2 >> window length)
            ]';  % 7 x nSeeds

        model_fit = @pid_biphasic_out;
        options = optimoptions('lsqnonlin', ...
            'MaxFunctionEvaluations',1e5, 'MaxIterations',1e5, ...
            'Display','off');

        paramsout = zeros(nParams, nOdors);
        resnormPerOdor = zeros(1, nOdors);
        exitflag = zeros(1, nOdors);
        output = cell(1, nOdors);
        residual = cell(1, nOdors);

        for k = 1:nOdors
            nBadX = sum(isnan(xdata(:,k)) | isinf(xdata(:,k)));
            nBadY = sum(isnan(ydata(:,k)) | isinf(ydata(:,k)));
            if nBadX > 0 || nBadY > 0
                warning('Odor %d: %d bad xdata samples, %d bad ydata samples -- skipping fit for this odor.', k, nBadX, nBadY);
                continue
            end

            bestResnorm = Inf;
            for s = 1:size(seedList,2)
                P0 = seedList(:,s);
                P0 = min(max(P0, lb1), ub1); % keep seed within bounds

                try
                    [Pfit, rn, res, ef, out] = lsqnonlin( ...
                        @(P) biphasic_residuals_single(P, k), P0, lb1, ub1, options);
                catch ME
                    warning('Odor %d, seed %d failed: %s', k, s, ME.message);
                    continue
                end

                if rn < bestResnorm
                    bestResnorm       = rn;
                    paramsout(:,k)    = Pfit;
                    resnormPerOdor(k) = rn;
                    exitflag(k)       = ef;
                    output{k}         = out;
                    residual{k}       = res;
                end
            end
        end
        resnorm = sum(resnormPerOdor);

        [kernelsout,PIDOut] = makeKenel(ValveVector,paramsout,timestep);

end


%% 2 : define the convolution function (model_fit)
    function [zdata] = pid_out(StartingKernels,xdata)

        % make loal copies
        Starting_Kernels    = StartingKernels;
        x_data              = xdata;
        maxT                = size(x_data,1)+1;
        convmode = 'full'; %'same';
        for kk = 1:size(x_data,2) 
            zdata(:,kk) = conv(x_data(:,kk)',Starting_Kernels(:,kk)',convmode);
        end

        zdata(maxT:end,:) = [];
    end

%% 2: parametric kernel + convolution
function [zdata] = pid_param_out(P, xdata)
    x_data = xdata;
    N      = size(x_data,1);
    dt     = timestep;  % capture from outer scope (nested function)
    tvec   = (0:N-1)' * dt;  % kernel time axis, same length as data

    zdata = zeros(N, size(x_data,2));
    for kk = 1:size(x_data,2)
        A     = P(1,kk);
        tRise = P(2,kk);
        tDecay= P(3,kk);
        delay = P(4,kk);

        t = tvec - delay;
        h = A * (exp(-t/tDecay) - exp(-t/tRise));
        h(t < 0) = 0;   % causal kernel

        full = conv(x_data(:,kk)', h', 'full');
        zdata(:,kk) = full(1:N)';
    end
end

function [zdata] = pid_adapt_out(P, xdata)
    x_data = xdata;
    N      = size(x_data,1);
    dt     = timestep;              % captured from outer scope
    tvec   = (0:N-1)' * dt;

    zdata = zeros(N, size(x_data,2));
    for kk = 1:size(x_data,2)
        A0    = P(1,kk);
        Ass   = P(2,kk);
        tauR  = P(3,kk);
        tauA  = P(4,kk);
        delay = P(5,kk);

        t = tvec - delay;
        S = (Ass + (A0-Ass).*exp(-t./tauA)) .* (1 - exp(-t./tauR));
        S(t < 0) = 0;                % causal: no response before delay

        h = [0; diff(S)] / dt;       % kernel = derivative of step response

        full = conv(x_data(:,kk)', h', 'full');
        zdata(:,kk) = full(1:N)';
    end
end

function [res] = biphasic_residuals_single(P, k)
    % Residual vector for lsqnonlin, single odor k: data-fit error
    % stacked with a regularization penalty on A2 (element 4 of P),
    % which suppresses the negative/adapting lobe unless it
    % meaningfully improves the fit for this odor.
    zdata = pid_biphasic_out(P, xdata(:,k));
    dataResidual = zdata - ydata(:,k);

    penalty = lambda * P(4);
    res = [dataResidual(:); penalty];
end

function [zdata] = pid_biphasic_out(P, xdata)
    x_data = xdata;
    N      = size(x_data,1);
    dt     = timestep;              % captured from outer scope
    tvec   = (0:N-1)' * dt;

    zdata = zeros(N, size(x_data,2));
    for kk = 1:size(x_data,2)
        A1    = P(1,kk); tau1 = P(2,kk); n1 = P(3,kk);
        A2    = P(4,kk); dtau = P(5,kk); n2 = P(6,kk);
        delay = P(7,kk);
        tau2  = tau1 + dtau; % ensures tau2 > tau1 without an absolute floor on tau2

        t = tvec - delay;
        t(t<=0) = NaN;   % guard against 0^n / log issues, zero out after

        g1 = (t/tau1).^n1 .* exp(n1 - n1*t/tau1);
        g2 = (t/tau2).^n2 .* exp(n2 - n2*t/tau2);
        g1(isnan(g1)) = 0;
        g2(isnan(g2)) = 0;

        h = A1*g1 - A2*g2;

        full = conv(x_data(:,kk)', h', 'full');
        zdata(:,kk) = full(1:N)';
    end
end
%%
    function [K,PIDOut] = makeKenel(ValveVector,Params,timestep,N)
        if nargin<4
            N = 5/timestep;
        end
        tvec   = (0:N-1)' * timestep;
        for kk = 1:size(Params,2)
            switch size(Params,1)
                case 4
                    A     = Params(1,kk);
                    tRise = Params(2,kk);
                    tDecay= Params(3,kk);
                    delay = Params(4,kk);
                    t = tvec - delay;
                    h = A * (exp(-t/tDecay) - exp(-t/tRise));
                    h(t < 0) = 0;   % causal kernel
                    K(:,kk) = h;
                case 5
                    A0    = Params(1,kk);
                    Ass   = Params(2,kk);
                    tauR  = Params(3,kk);
                    tauA  = Params(4,kk);
                    delay = Params(5,kk);
                    t = tvec - delay;
                    S = (Ass + (A0-Ass).*exp(-t./tauA)) .* (1 - exp(-t./tauR));
                    S(t < 0) = 0;                % causal: no response before delay
                    h = [0; diff(S)] / timestep;       % kernel = derivative of step response
                    K(:,kk) = h;
                case 7
                    A1    = Params(1,kk); 
                    tau1 = Params(2,kk); 
                    n1 = Params(3,kk);
                    A2    = Params(4,kk); 
                    dtau = Params(5,kk); 
                    n2 = Params(6,kk);
                    delay = Params(7,kk);
                    tau2 = tau1 + dtau; % ensures tau2 > tau1 without an absolute floor on tau2
                    t = tvec - delay;
                    t(t<=0) = NaN;   % guard against 0^n / log issues, zero out after
                    g1 = (t/tau1).^n1 .* exp(n1 - n1*t/tau1);
                    g2 = (t/tau2).^n2 .* exp(n2 - n2*t/tau2);
                    g1(isnan(g1)) = 0;
                    g2(isnan(g2)) = 0;
                    h = A1*g1 - A2*g2;
                    K(:,kk) = h;
            end

            full = conv(ValveVector,h','full');
            PIDOut(:,kk) = full(1:numel(ValveVector))';
        end
    end


%% stop clock
toc

end