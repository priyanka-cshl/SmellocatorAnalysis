function Fitted = FitCIDPIDKernels_Variants(PIDdir, varargin)
% Fits a kernel that, convolved with the odor valve waveform, recreates
% the measured PID waveform for each odor.
%
% Model: PID(t) = kernel(t) * drive(t), where:
%   - drive(t) = valve(t) .* [Ainf + (1-Ainf)*exp(-t/tau_source)]
%     A signed source-dynamics term: the odor source's availability
%     relaxes from 1 (at valve-open) toward Ainf. Ainf<1 = depletion
%     (headspace runs low during the pulse), Ainf>1 = buildup (source
%     takes time to reach full strength), Ainf=1 = no source dynamics.
%   - kernel(t) = A * (exp(-(t/tau_decay)^beta) - exp(-t/tau_rise))
%     A stretched-exponential difference-of-exponentials kernel:
%     tau_rise governs the fast onset, tau_decay/beta govern the
%     (possibly non-exponential/fat-tailed) approach to plateau and
%     offset decay. beta=1 is a plain exponential tail; beta<1 gives a
%     fatter, slower-converging tail. tau_decay is internally
%     parametrized as tau_rise + extra_decay (extra_decay >= a small
%     floor) to guarantee tau_decay > tau_rise and avoid the
%     near-cancellation degeneracy that occurs when the two time
%     constants collapse together.
%
% Fit via multi-start lsqnonlin per odor (7 free parameters), with a
% light shrinkage penalty pulling beta toward 1 -- this keeps
% low-SNR/noisy odors from using the extra beta flexibility to
% overfit noise (which shows up as spurious ringing at the onset/
% offset transitions); beta only deviates from 1 when the data
% clearly supports it.
%
% Inputs: Folder with the processed PID data.
% from averagedPID.mat: PIDTraces (n x t), n = odors, t = samples @500hz,
% and window - time before and after odor off.
% from quickprocessPID.mat: stimulus timing info, used to build the
% valve waveform.

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 10, @(x) isnumeric(x)); % in ms

params.parse(varargin{:});
binsize = params.Results.binsize;

PID     = load(fullfile(PIDdir,'quickprocessPID.mat'));
Traces  = load(fullfile(PIDdir,'averagedPID.mat'));

OdorDuration = PID.StimSettings.timing(3)/1000;
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

tic

%% Fit: signed source-dynamics drive x stretched-exponential kernel
% params per odor: [A, tau_rise, extra_decay, beta, Ainf, tau_source, delay]
% (tau_decay = tau_rise + extra_decay internally)
nParams = 7;

lb1 = [0,    0.005, 0.02, 0.2, 0,   0.05, 0]';
ub1 = [5,    1,     3,    1.0, 3,   20,   0.5]';

seedList = [ ...
    1,    0.05, 0.30, 1.0, 1.0, 20,   0;      % fast onset, no source dynamics, plain exponential
    1,    0.05, 0.30, 0.5, 0.3, 0.5,  0;      % fast onset, depletion, moderate stretch
    1,    0.05, 0.30, 0.5, 1.0, 20,   0;      % fast onset, no source dynamics, stretched
    1,    0.10, 0.50, 0.5, 2.0, 1.0,  0;      % buildup, stretched
    1,    0.30, 0.80, 0.6, 1.0, 20,   0;      % slow onset (kernel-driven), no source dynamics
    0.10, 0.01, 0.17, 1.0, 0.57, 0.05, 0.16;  % fast jump + fast depletion, plain exponential
    0.10, 0.01, 0.30, 0.5, 0.57, 0.05, 0.16;  % fast jump + fast depletion, stretched
    ]';

lambdaBeta = 0.05; % shrinkage weight pulling beta toward 1 -- larger = more conservative (less prone to overfitting noise)

options = optimoptions('lsqnonlin', ...
    'MaxFunctionEvaluations',1e5, 'MaxIterations',1e5, ...
    'Display','off');

paramsout = zeros(nParams, nOdors);
resnormPerOdor = zeros(1, nOdors);
exitflag = zeros(1, nOdors);
output = cell(1, nOdors);
residual = cell(1, nOdors);
kernelsout = zeros(size(xdata,1), nOdors);
PIDOut     = zeros(size(xdata,1), nOdors);
DriveOut   = zeros(size(xdata,1), nOdors);

for k = 1:nOdors
    nBadX = sum(isnan(xdata(:,k)) | isinf(xdata(:,k)));
    nBadY = sum(isnan(ydata(:,k)) | isinf(ydata(:,k)));
    if nBadX > 0 || nBadY > 0
        warning('Odor %d: %d bad xdata samples, %d bad ydata samples -- skipping fit for this odor.', k, nBadX, nBadY);
        continue
    end

    bestResnorm = Inf;
    for s = 1:size(seedList,2)
        P0 = min(max(seedList(:,s), lb1), ub1);

        try
            [Pfit, rn, res, ef, out] = lsqnonlin( ...
                @(P) stretched_residuals_single(P, k), P0, lb1, ub1, options);
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

    [K,PO,DO] = makeStretchedKernel(ValveVector, paramsout(:,k), timestep, TimeVector);
    kernelsout(:,k) = K; PIDOut(:,k) = PO; DriveOut(:,k) = DO;
end
resnorm = sum(resnormPerOdor);

fprintf('\n--- Fit summary (A, tau_rise, tau_decay, beta, Ainf, tau_source, delay) ---\n');
for k = 1:nOdors
    if all(paramsout(:,k)==0), continue; end
    tau_decay = paramsout(2,k) + paramsout(3,k);
    fprintf('Odor %2d: A=%.3f  tau_rise=%.4f  tau_decay=%.4f  beta=%.3f  Ainf=%.3f  tau_source=%.3f  delay=%.3f\n', ...
        k, paramsout(1,k), paramsout(2,k), tau_decay, paramsout(4,k), paramsout(5,k), paramsout(6,k), paramsout(7,k));
end
fprintf('--- end summary ---\n\n');

Fitted = struct('kernelsout', kernelsout, 'PIDOut', PIDOut, 'DriveOut', DriveOut, 'paramsout', paramsout, ...
    'TimeVector', TimeVector, 'ValveVector', ValveVector, 'PIDtraces', PIDtraces);

%% Model functions

    function [res] = stretched_residuals_single(P, k)
        zdata = pid_stretched_out(P, xdata(:,k));
        dataResidual = zdata - ydata(:,k);
        betaPenalty = lambdaBeta * (1 - P(4)); % shrink beta toward 1 unless data justifies deviating
        res = [dataResidual(:); betaPenalty];
    end

    function [zdata] = pid_stretched_out(P, xdata_col)
        A           = P(1);
        tau_rise    = P(2);
        extra_decay = P(3);
        beta        = P(4);
        Ainf        = P(5);
        tau_source  = P(6);
        delay       = P(7);
        tau_decay   = tau_rise + extra_decay; % ensures tau_decay > tau_rise, avoids near-cancellation degeneracy

        N  = size(xdata_col,1);
        dt = timestep;
        tvec = (0:N-1)' * dt;

        ageSinceOnset  = TimeVector(:);
        sourceFactor   = Ainf + (1-Ainf) * exp(-max(ageSinceOnset,0)/tau_source);
        effectiveDrive = xdata_col .* sourceFactor;

        t = tvec - delay;
        t(t < 0) = NaN;  % guard against negative base ^ non-integer beta
        h = A * (exp(-(t/tau_decay).^beta) - exp(-t/tau_rise));
        h(isnan(h)) = 0;

        full = conv(effectiveDrive', h', 'full');
        zdata = full(1:N)';
    end

    function [K,PIDOut_,DriveOut_] = makeStretchedKernel(ValveVector_,P,timestep_,TimeVectorIn)
        N  = numel(ValveVector_);
        dt = timestep_;
        tvec = (0:N-1)' * dt;
        ageSinceOnset = TimeVectorIn(:);

        A           = P(1);
        tau_rise    = P(2);
        extra_decay = P(3);
        beta        = P(4);
        Ainf        = P(5);
        tau_source  = P(6);
        delay       = P(7);
        tau_decay   = tau_rise + extra_decay;

        sourceFactor = Ainf + (1-Ainf) * exp(-max(ageSinceOnset,0)/tau_source);
        DriveOut_ = ValveVector_(:) .* sourceFactor;

        t = tvec - delay;
        t(t < 0) = NaN;
        K = A * (exp(-(t/tau_decay).^beta) - exp(-t/tau_rise));
        K(isnan(K)) = 0;

        full = conv(DriveOut_', K', 'full');
        PIDOut_ = full(1:N)';
    end

%% Joint fit comparison: test which parameters can be shared across
% odors without hurting the fit. Confirmed already: delay, tau_rise,
% tau_decay share with ~0% cost (rig/plumbing properties). Here we
% additionally test forcing Ainf or tau_source to be shared too, to see
% whether source-depletion behavior is also a shared rig property or
% genuinely odor-specific (chemistry-dependent).
goodOdors = find(any(paramsout~=0,1));  % exclude any odors that failed to fit at all
paramNames = {'A','tau_rise','extra_decay','beta','Ainf','tau_source','delay'};

baseSharedIdx = [7 2 3];  % delay, tau_rise, extra_decay (already confirmed ~0% cost)

configs = { ...
    'baseline (delay+rise+decay shared)',        baseSharedIdx; ...
    'baseline + shared Ainf',                     [baseSharedIdx 5]; ...
    'baseline + shared tau_source',                [baseSharedIdx 6]; ...
    };

optionsJoint = optimoptions('lsqnonlin', ...
    'MaxFunctionEvaluations',1e6, 'MaxIterations',1e6, ...
    'Display','off');

fprintf('\n--- Joint-fit configuration comparison ---\n');
fprintf('Independent per-odor total resnorm: %.4f\n', resnorm);

Fitted.jointConfigs = struct('name',{},'sharedIdx',{},'sharedFit',{},'resnorm',{}, ...
    'pctCost',{},'paramsoutJoint',{},'kernelsoutJoint',{},'PIDOutJoint',{},'DriveOutJoint',{});

for c = 1:size(configs,1)
    cName = configs{c,1};
    sIdx  = configs{c,2};
    pIdx  = setdiff(1:7, sIdx);

    sInit = median(paramsout(sIdx, goodOdors), 2);
    pInit = paramsout(pIdx, goodOdors);
    P0c = [sInit; pInit(:)];

    lbSh = lb1(sIdx); ubSh = ub1(sIdx);
    lbPo = repmat(lb1(pIdx),1,numel(goodOdors));
    ubPo = repmat(ub1(pIdx),1,numel(goodOdors));
    lbC = [lbSh; lbPo(:)];
    ubC = [ubSh; ubPo(:)];

    [PfitC, rnC] = lsqnonlin( ...
        @(Pj) joint_residuals_generic(Pj, sIdx, pIdx, goodOdors), ...
        P0c, lbC, ubC, optionsJoint);

    nShC = numel(sIdx);
    sharedFitC  = PfitC(1:nShC);
    perOdorFitC = reshape(PfitC(nShC+1:end), numel(pIdx), numel(goodOdors));

    paramsoutC  = zeros(7, nOdors);
    kernelsoutC = zeros(size(xdata,1), nOdors);
    PIDOutC     = zeros(size(xdata,1), nOdors);
    DriveOutC   = zeros(size(xdata,1), nOdors);
    for kk = 1:numel(goodOdors)
        Pk = zeros(7,1);
        Pk(sIdx) = sharedFitC;
        Pk(pIdx) = perOdorFitC(:,kk);
        paramsoutC(:,goodOdors(kk)) = Pk;

        [Kj,POj,DOj] = makeStretchedKernel(ValveVector, Pk, timestep, TimeVector);
        kernelsoutC(:,goodOdors(kk)) = Kj;
        PIDOutC(:,goodOdors(kk))     = POj;
        DriveOutC(:,goodOdors(kk))   = DOj;
    end

    pctCost = 100*(rnC-resnorm)/resnorm;

    % Per-odor cost: the pooled total above is dominated by the 3 noisy
    % odors (their individual resnorm is ~100-1000x the clean odors'),
    % so a small pooled % change could still hide a large relative cost
    % for individual well-behaved odors. Compute per-odor resnorm under
    % this config and compare directly to that odor's own independent
    % fit resnorm.
    resnormPerOdorC = nan(1,nOdors);
    for kk = 1:numel(goodOdors)
        oi = goodOdors(kk);
        resnormPerOdorC(oi) = sum((PIDOutC(:,oi) - ydata(:,oi)).^2);
    end
    pctCostPerOdor = 100*(resnormPerOdorC - resnormPerOdor) ./ resnormPerOdor;

    fprintf('[%s]: shared={%s} -> resnorm=%.4f (%.1f%% pooled cost)\n', ...
        cName, strjoin(paramNames(sIdx),','), rnC, pctCost);
    [worstCost, worstOdor] = max(pctCostPerOdor(goodOdors));
    fprintf('    worst-affected odor: #%d, %.1f%% increase in that odor''s own resnorm\n', ...
        goodOdors(worstOdor), worstCost);
    fprintf('    median per-odor cost: %.1f%%\n', median(pctCostPerOdor(goodOdors)));

    Fitted.jointConfigs(c).name            = cName;
    Fitted.jointConfigs(c).sharedIdx       = sIdx;
    Fitted.jointConfigs(c).sharedFit       = sharedFitC;
    Fitted.jointConfigs(c).resnorm         = rnC;
    Fitted.jointConfigs(c).pctCost         = pctCost;
    Fitted.jointConfigs(c).resnormPerOdor  = resnormPerOdorC;
    Fitted.jointConfigs(c).pctCostPerOdor  = pctCostPerOdor;
    Fitted.jointConfigs(c).paramsoutJoint  = paramsoutC;
    Fitted.jointConfigs(c).kernelsoutJoint = kernelsoutC;
    Fitted.jointConfigs(c).PIDOutJoint     = PIDOutC;
    Fitted.jointConfigs(c).DriveOutJoint   = DriveOutC;
end
fprintf('--- end comparison ---\n\n');

    function res = joint_residuals_generic(Pj, sIdxArg, pIdxArg, goodOdorsArg)
        % All varying inputs are explicit function arguments (not
        % implicitly shared with the parent workspace), so this is safe
        % to reuse across multiple configurations with no risk of
        % cross-call variable collisions.
        nShG = numel(sIdxArg);
        shG = Pj(1:nShG);
        poG = reshape(Pj(nShG+1:end), numel(pIdxArg), numel(goodOdorsArg));
        res = [];
        for kkk = 1:numel(goodOdorsArg)
            oidx_g  = goodOdorsArg(kkk);
            Pfull_g = zeros(7,1);
            Pfull_g(sIdxArg) = shG;
            Pfull_g(pIdxArg) = poG(:,kkk);
            zdata_g = pid_stretched_out(Pfull_g, xdata(:,oidx_g));
            dataResidual_g = zdata_g - ydata(:,oidx_g);
            betaPenalty_g = lambdaBeta * (1 - Pfull_g(4));
            res = [res; dataResidual_g(:); betaPenalty_g]; %#ok<AGROW>
        end
    end

toc

end