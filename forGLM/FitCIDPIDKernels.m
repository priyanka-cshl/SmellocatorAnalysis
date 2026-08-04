function Fitted = FitCIDPIDKernels(PIDdir, varargin)
% Fits a two-regime model that, given the odor valve waveform,
% recreates the measured PID waveform for each odor.
%
% ON regime (valve open):
%   PID(t) = kernel_on(t) * drive(t)
%   drive(t)     = valve(t) .* [Ainf + (1-Ainf)*exp(-t/tau_source)]
%                  (source-depletion: odor availability relaxes toward Ainf)
%   kernel_on(t) = A * [exp(-((t-delay)/tau_decay_on)^beta_on) - exp(-(t-delay)/tau_rise)]
%                  (stretched difference-of-exponentials: onset speed + approach to plateau)
%
% OFF regime (valve closed): two independent clearance processes,
% anchored to wherever the ON regime left off at the moment the valve
% closure is felt at the sensor (shifted by delay_off, a separate lag
% from the onset delay):
%   PID(t) = PID(t_off) * [(1-w_slow)*exp(-(t-t_off)/tau_rise) + w_slow*exp(-(t-t_off)/tau_decay_off_slow)]
%   - fast fraction (1-w_slow) reuses tau_rise (mirrors the onset transport)
%   - slow fraction w_slow has its own tau_decay_off_slow (odor lingering/diffusing out of tubing)
%
% Design notes (why it looks like this):
%   - ON and OFF are governed by separate parameters because a single
%     kernel forcing one timescale to fit both "rise to plateau" and
%     "decay after close" was structurally unable to do either well.
%   - The OFF decay is a two-component mixture (fast + slow) rather
%     than one curve because the data looks like the sum of two
%     physical processes (fast transport clearing + slow tubing
%     residue), not one non-exponential shape.
%   - delay, delay_off, and w_slow are shared across all odors.
%     delay/delay_off are rig/plumbing properties (confirmed: merging
%     them into one lag costs ~99% median per-odor fit quality --
%     valve open/close timing genuinely differs). w_slow (the fast/slow
%     OFF-clearance mixing fraction) was tested empirically and found
%     shareable at negligible cost (~0%), suggesting it's a rig/
%     geometry property rather than odor chemistry. tau_rise is
%     deliberately NOT shared even though it's reused in both regimes:
%     sharing it forced one global value to trade off ON-fit quality
%     against OFF-fit quality across all 16 odors at once. Letting it
%     vary per odor removed that compromise.
%   - Everything else (A, tau_decay_on, beta_on, Ainf, tau_source,
%     tau_decay_off_slow) is fit per odor -- tau_decay_off_slow was
%     tested for sharing too and found genuinely odor-specific
%     (~91% worst-case cost), consistent with it reflecting odor
%     stickiness/surface chemistry rather than rig geometry.
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

for i = 1:nOdors
    PIDfiltered(:,i) = PIDfiltered(:,i) - mean(PIDfiltered(1:500*abs(window(1)),i));
    PIDfiltered(:,i) = PIDfiltered(:,i)/max(PIDfiltered(:,i));
end

[PIDtraces] = interp1(Traces.TimeOut(1:sampMax,1),PIDfiltered,TimeVector);
PIDtraces(1:startat,:) = [];
TimeVector(:,1:startat) = [];
ValveVector(:,1:startat) = [];

nOdors = size(PIDtraces,2);
xdata = repmat(ValveVector',1,nOdors);
ydata = PIDtraces;

%% Step 1: independent per-odor fit (multi-start), used as warm start
% for the shared-parameter fit below.
% params per odor: [A, tau_rise, tau_decay_on, beta_on, Ainf, tau_source, delay, tau_decay_off_slow, w_slow, delay_off]
paramNames = {'A','tau_rise','tau_decay_on','beta_on','Ainf','tau_source','delay','tau_decay_off_slow','w_slow','delay_off'};

lb1 = [0,    0.005, 0.02, 0.2, 0,   0.05, 0,   0.05, 0, 0  ]';
ub1 = [5,    1,     3,    1.0, 3,   20,   0.5, 10,   1, 0.5]';

seedList = [ ...
    1,    0.05, 0.30, 1.0, 1.0, 20,   0,    1.0, 0.3, 0.15;
    1,    0.05, 0.30, 0.5, 0.3, 0.5,  0,    2.0, 0.3, 0.15;
    1,    0.05, 0.30, 0.5, 1.0, 20,   0,    0.5, 0.5, 0.15;
    1,    0.10, 0.50, 0.5, 2.0, 1.0,  0,    1.5, 0.4, 0.15;
    1,    0.30, 0.80, 0.6, 1.0, 20,   0,    2.0, 0.5, 0.15;
    0.10, 0.01, 0.17, 1.0, 0.57, 0.05, 0.16, 1.0, 0.3, 0.15;
    0.10, 0.01, 0.30, 0.5, 0.57, 0.05, 0.16, 2.0, 0.4, 0.15;
    1,    0.05, 0.30, 1.0, 1.0, 20,   0,    3.0, 0.5, 0.15;   % strongly dominated by slow tubing-clearance
    1,    0.05, 0.30, 0.5, 1.0, 20,   0,    0.3, 0.1, 0.15;   % mostly fast-clearing, small slow tail
    1,    0.05, 0.30, 1.0, 1.0, 20,   0,    1.0, 0.3, 0.30;   % offset delay noticeably longer than onset delay
    ]';

lambdaBeta = 0.05; % shrinkage weight pulling beta_on toward 1 -- keeps noisy odors from overfitting

options = optimoptions('lsqnonlin', ...
    'MaxFunctionEvaluations',1e5, 'MaxIterations',1e5, ...
    'Display','off');

nParams = 10;
paramsoutIndep = zeros(nParams, nOdors);
for k = 1:nOdors
    nBadX = sum(isnan(xdata(:,k)) | isinf(xdata(:,k)));
    nBadY = sum(isnan(ydata(:,k)) | isinf(ydata(:,k)));
    if nBadX > 0 || nBadY > 0
        warning('Odor %d: %d bad xdata samples, %d bad ydata samples -- skipping.', k, nBadX, nBadY);
        continue
    end

    bestResnorm = Inf;
    for s = 1:size(seedList,2)
        P0 = min(max(seedList(:,s), lb1), ub1);
        try
            [Pfit, rn] = lsqnonlin(@(P) tworegime_residuals_single(P, k), P0, lb1, ub1, options);
        catch
            continue
        end
        if rn < bestResnorm
            bestResnorm = rn;
            paramsoutIndep(:,k) = Pfit;
        end
    end
end

%% Step 2: joint fit with delay and delay_off shared across all odors
% (tau_rise deliberately left per-odor -- see header notes), using
% Step 1 as warm start.
goodOdors = find(any(paramsoutIndep~=0,1));
sharedIdx  = [7 10 9];  % delay, delay_off, w_slow
perOdorIdx = setdiff(1:nParams, sharedIdx);

sharedInit  = median(paramsoutIndep(sharedIdx, goodOdors), 2);
perOdorInit = paramsoutIndep(perOdorIdx, goodOdors);
Pjoint0 = [sharedInit; perOdorInit(:)];

nShared = numel(sharedIdx);
lbPerOdorRep = repmat(lb1(perOdorIdx), 1, numel(goodOdors));
ubPerOdorRep = repmat(ub1(perOdorIdx), 1, numel(goodOdors));
lbJoint = [lb1(sharedIdx); lbPerOdorRep(:)];
ubJoint = [ub1(sharedIdx); ubPerOdorRep(:)];

optionsJoint = optimoptions('lsqnonlin', ...
    'MaxFunctionEvaluations',1e6, 'MaxIterations',1e6, ...
    'Display','off');

[PjointFit, resnormJoint] = lsqnonlin( ...
    @(Pj) joint_residuals(Pj, sharedIdx, perOdorIdx, goodOdors), ...
    Pjoint0, lbJoint, ubJoint, optionsJoint);

sharedFit  = PjointFit(1:nShared);
perOdorFit = reshape(PjointFit(nShared+1:end), numel(perOdorIdx), numel(goodOdors));

paramsout   = zeros(nParams, nOdors);
kernelsout  = zeros(size(xdata,1), nOdors);
PIDOut      = zeros(size(xdata,1), nOdors);
DriveOut    = zeros(size(xdata,1), nOdors);
for kk = 1:numel(goodOdors)
    Pk = zeros(nParams,1);
    Pk(sharedIdx)  = sharedFit;
    Pk(perOdorIdx) = perOdorFit(:,kk);
    paramsout(:,goodOdors(kk)) = Pk;

    [K,PO,DO] = make_tworegime_output(Pk, xdata(:,goodOdors(kk)));
    kernelsout(:,goodOdors(kk)) = K;
    PIDOut(:,goodOdors(kk))     = PO;
    DriveOut(:,goodOdors(kk))   = DO;
end

fprintf('\nShared parameters (fit jointly across all odors):\n');
fprintf('  delay      = %.4f s\n', sharedFit(1));
fprintf('  delay_off  = %.4f s\n', sharedFit(2));
fprintf('  w_slow     = %.4f\n', sharedFit(3));
fprintf('  (tau_rise, tau_decay_off_slow, and all other params are free per-odor)\n');
fprintf('Total resnorm: %.4f\n\n', resnormJoint);

%% Sanity check
resnormPerOdor = nan(1,nOdors);
boundHits = {};
for k = goodOdors
    resnormPerOdor(k) = sum((PIDOut(:,k) - ydata(:,k)).^2);
    for p = 1:nParams
        tol = 0.005*(ub1(p)-lb1(p));
        if abs(paramsout(p,k)-lb1(p)) < tol || abs(paramsout(p,k)-ub1(p)) < tol
            boundHits{end+1} = sprintf('odor %d: %s pinned at bound (%.4g)', k, paramNames{p}, paramsout(p,k)); %#ok<AGROW>
        end
    end
end

fprintf('--- Sanity check ---\n');
fprintf('Per-odor resnorm: median = %.4f, max = %.4f (odor %d)\n', ...
    median(resnormPerOdor(goodOdors)), max(resnormPerOdor(goodOdors)), ...
    goodOdors(resnormPerOdor(goodOdors)==max(resnormPerOdor(goodOdors))));
if isempty(boundHits)
    fprintf('No parameters pinned at their bounds.\n');
else
    fprintf('%d parameter(s) pinned at bounds (worth checking if that bound is too tight):\n', numel(boundHits));
    for i = 1:numel(boundHits)
        fprintf('  %s\n', boundHits{i});
    end
end
if any(PIDOut(:) > 1.3) || any(PIDOut(:) < -0.3)
    warning('Model predicts values outside [-0.3, 1.3] somewhere -- check for overshoot/undershoot artifacts.');
else
    fprintf('No overshoot/undershoot artifacts detected (predictions stay within [-0.3, 1.3]).\n');
end
fprintf('--- end sanity check ---\n\n');

% Parameter-sharing questions already settled empirically (per-odor
% cost testing, not just pooled resnorm):
%   - delay, delay_off: must stay separate (merging costs ~99% median,
%     ~300% worst-case per-odor) -- valve open/close timing differs.
%   - w_slow: safe to share (~0-2% cost) -- adopted above.
%   - tau_decay_off_slow: NOT safe to share (~91% worst-case cost) --
%     genuinely odor-specific, kept per-odor.
%   - tau_rise: NOT safe to share -- reused in both ON and OFF regimes,
%     sharing it forces one value to trade off ON-fit vs OFF-fit
%     quality across all odors at once (confirmed via direct
%     before/after comparison). Kept per-odor.

Fitted = struct('kernelsout', kernelsout, 'PIDOut', PIDOut, 'DriveOut', DriveOut, ...
    'paramsout', paramsout, 'sharedParams', sharedFit, 'resnormPerOdor', resnormPerOdor, ...
    'TimeVector', TimeVector, 'ValveVector', ValveVector, 'PIDtraces', PIDtraces);

%% Model functions

    function [res] = tworegime_residuals_single(P, k)
        zdata = pid_tworegime_out(P, xdata(:,k));
        dataResidual = zdata - ydata(:,k);
        betaPenalty = lambdaBeta * (1 - P(4)); % shrink beta_on toward 1
        res = [dataResidual(:); betaPenalty];
    end

    function res = joint_residuals(Pj, sIdxArg, pIdxArg, goodOdorsArg)
        nShG = numel(sIdxArg);
        shG = Pj(1:nShG);
        poG = reshape(Pj(nShG+1:end), numel(pIdxArg), numel(goodOdorsArg));
        res = [];
        for kkk = 1:numel(goodOdorsArg)
            oidx_g  = goodOdorsArg(kkk);
            Pfull_g = zeros(nParams,1);
            Pfull_g(sIdxArg) = shG;
            Pfull_g(pIdxArg) = poG(:,kkk);
            zdata_g = pid_tworegime_out(Pfull_g, xdata(:,oidx_g));
            dataResidual_g = zdata_g - ydata(:,oidx_g);
            betaPenalty_g = lambdaBeta * (1 - Pfull_g(4));
            res = [res; dataResidual_g(:); betaPenalty_g]; %#ok<AGROW>
        end
    end


    function [zdata] = pid_tworegime_out(P, xdata_col)
        [zdata, ~, ~] = compute_tworegime(P, xdata_col);
    end

    function [K,PIDOut_,DriveOut_] = make_tworegime_output(P, xdata_col)
        [PIDOut_, K, DriveOut_] = compute_tworegime(P, xdata_col);
    end

    function [zdata, K, DriveOut_] = compute_tworegime(P, xdata_col)
        A            = P(1);
        tau_rise     = P(2);
        tau_decay_on = P(3);
        beta_on      = P(4);
        Ainf         = P(5);
        tau_source   = P(6);
        delay        = P(7);
        tau_decay_off_slow = P(8);
        w_slow       = P(9);
        delay_off    = P(10);

        N  = size(xdata_col,1);
        dt = timestep;
        tvec = (0:N-1)' * dt;

        ageSinceOnset  = TimeVector(:);
        sourceFactor   = Ainf + (1-Ainf) * exp(-max(ageSinceOnset,0)/tau_source);
        DriveOut_ = xdata_col .* sourceFactor;

        % ON regime: same drive-convolution kernel as before
        t = tvec - delay;
        t(t < 0) = NaN;
        K = A * (exp(-(t/tau_decay_on).^beta_on) - exp(-t/tau_rise));
        K(isnan(K)) = 0;

        full = conv(DriveOut_', K', 'full');
        y_on = full(1:N)';

        % OFF regime: TWO independent clearance processes, anchored to
        % the ON regime's value at the moment the valve closure is
        % actually felt at the sensor -- shifted by delay_off, a
        % SEPARATE transport lag from the onset delay (valve open/close
        % mechanics can differ). A fraction (1-w_slow) clears fast,
        % reusing this odor's own tau_rise; a fraction w_slow clears
        % slowly via its own independent tau_decay_off_slow (odor
        % lingering/diffusing out of tubing).
        rawSwitchIdx = find(xdata_col==1, 1, 'last');
        zdata = y_on;
        if ~isempty(rawSwitchIdx)
            switchIdx = rawSwitchIdx + round(delay_off/dt);
            switchIdx = min(switchIdx, N);
            if switchIdx < N
                y_boundary = y_on(switchIdx);
                tOff = tvec(switchIdx+1:end) - tvec(switchIdx);
                fastPart = (1-w_slow) * exp(-tOff/tau_rise);
                slowPart = w_slow     * exp(-tOff/tau_decay_off_slow);
                zdata(switchIdx+1:end) = y_boundary * (fastPart + slowPart);
            end
        end
    end

end