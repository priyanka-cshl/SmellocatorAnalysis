function [Fitted,resnorm,residual,exitflag,output] = FitCIDPIDKernels3(PIDdir, varargin)
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
KernelMode = 'stretched_test';

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

    case 'basis'
        % Raised-cosine log-time basis (Pillow-style), fit via linear
        % least squares with a SMOOTHNESS penalty (not plain ridge).
        % Plain L2-on-weights (ridge) shrinks everything uniformly,
        % which can (a) still allow Gibbs-ringing notches at sharp
        % onset/offset transitions since it doesn't penalize
        % wiggliness specifically, and (b) suppress legitimate slow,
        % weak structure along with genuine noise. A roughness penalty
        % on the kernel's 2nd derivative instead discourages jagged
        % reconstructions while leaving smooth slow trends alone.
        nBasis    = 14;    % more bumps than before -- extra resolution for the slow-decay range
        kernelDur = 4;     % seconds of kernel support
        nK = round(kernelDur/timestep);
        tK = (0:nK-1)'*timestep;

        dt_min   = 2*timestep;      % pushed slightly later than 1*timestep to avoid an overly narrow first bump
        dt_max   = kernelDur*0.9;
        nlOffset = 1e-3;
        nlin     = @(x) log(x + nlOffset);

        ctr1   = nlin(dt_min);
        ctrEnd = nlin(dt_max);
        dCtr   = (ctrEnd - ctr1) / (nBasis - 1);
        ctrs   = ctr1:dCtr:ctrEnd;

        raisedCosine = @(t,ctr) (cos(max(-pi, min(pi, (nlin(t)-ctr)*pi/(2*dCtr)))) + 1) / 2;

        Basis = zeros(nK, nBasis);
        for b = 1:nBasis
            Basis(:,b) = raisedCosine(tK, ctrs(b));
        end

        % Taper the last portion of every basis function smoothly to
        % zero. Without this, the outermost basis function's own decay
        % can still be substantially nonzero right at nK (since dt_max
        % is only 90% of kernelDur), producing an abrupt truncation
        % discontinuity in the reconstructed kernel at its far edge.
        taperFrac = 0.15;
        taperLen  = round(taperFrac*nK);
        taperWin  = (cos(linspace(0,pi,taperLen))+1)/2;  % smooth 1 -> 0
        taperCol  = ones(nK,1);
        taperCol(end-taperLen+1:end) = taperWin;
        Basis = Basis .* taperCol;

        % Second-difference operator on the kernel time axis, used to
        % build a roughness penalty in basis-weight space: penalizing
        % D2*(Basis*w) = (D2*Basis)*w discourages wiggly/jagged kernel
        % reconstructions without directly shrinking the weights
        % themselves (so genuine slow, low-amplitude structure isn't
        % penalized just for being small).
        %
        % IMPORTANT: this penalty is TIME-VARYING, not uniform. A
        % uniform penalty was fighting itself: strong enough to smooth
        % out slow-decay ringing also flattened the legitimate fast
        % onset, causing the optimizer to ring/overshoot right at the
        % onset/offset transitions to compensate (the notches). Instead:
        % leave early kernel lags (~genuine fast rise) almost
        % unconstrained, and only enforce smoothness at longer lags
        % (where the response should be a gentle, slow trend, not noise).
        D2 = diff(eye(nK), 2, 1);      % (nK-2) x nK second-difference matrix

        earlyCutoff_sec = 0.3;                          % lags below this are left nearly free
        earlyN = round(earlyCutoff_sec/timestep);
        roughWeight = ones(nK-2,1);
        roughWeight(1:min(earlyN,end)) = 0.02;           % light penalty near onset (allow sharpness)
        Rough = (D2*Basis)' * diag(roughWeight) * (D2*Basis);  % nBasis x nBasis roughness matrix

        % Normalize Rough's scale to match the natural scale of the
        % weight-space Gram matrix (Basis'*Basis). Without this, Rough's
        % magnitude is set by an unnormalized discrete second-difference
        % (no division by timestep^2), which can be wildly different
        % from the data term U'U -- making lambdaSmooth's effective
        % strength unpredictable (in this case, almost certainly too
        % strong, swamping the fit everywhere despite the time-varying
        % weight). After normalizing, lambdaSmooth is a unitless
        % multiplier in a sane ~0.01-2 range.
        Rough = Rough / trace(Rough) * trace(Basis'*Basis);

        lambdaRidge = 1e-3;   % small ridge for numerical stability only
        lambdaSmooth = 0.3;   % main regularizer for the LATE (post-earlyCutoff) part of the kernel -- now unitless/normalized

        N = size(xdata,1);

        weightsout = zeros(nBasis, nOdors);
        kernelsout = zeros(nK, nOdors);
        PIDOut     = zeros(N, nOdors);
        resnormPerOdor = zeros(1, nOdors);

        for k = 1:nOdors
            U = zeros(N, nBasis);
            for b = 1:nBasis
                full = conv(xdata(:,k)', Basis(:,b)', 'full');
                U(:,b) = full(1:N)';
            end

            w = (U'*U + lambdaRidge*eye(nBasis) + lambdaSmooth*Rough) \ (U'*ydata(:,k));

            weightsout(:,k) = w;
            kernelsout(:,k) = Basis * w;
            PIDOut(:,k)     = U * w;
            resnormPerOdor(k) = sum((U*w - ydata(:,k)).^2);
        end
        resnorm  = sum(resnormPerOdor);
        exitflag = ones(1, nOdors);
        residual = PIDOut - ydata;
        output   = struct('nBasis', nBasis, 'kernelDur', kernelDur, ...
                           'basisCenters_sec', exp(ctrs), ...
                           'lambdaRidge', lambdaRidge, 'lambdaSmooth', lambdaSmooth);

    case 'depletion'
        % Physically-motivated alternative: rather than putting the sag
        % into the KERNEL shape (which needed a biphasic/negative lobe
        % and was hard to fit robustly), put it into the DRIVE signal.
        % Model: the odor headspace at the source depletes while the
        % valve is open (simple first-order depletion, time constant
        % tau_deplete), and this depleting drive is convolved with a
        % simple, always-positive kernel (fast rise, slower decay --
        % same form as the original 'exponentials' case). This
        % naturally reproduces fast rise -> plateau or gradual sag (if
        % tau_deplete is short relative to the pulse) -> offset decay,
        % with just 5 interpretable, well-behaved parameters per odor
        % and no need for negative lobes or basis functions.
        %
        % params per odor: [A, tau_rise, extra_decay, tau_deplete, delay]
        % note: tau_decay = tau_rise + extra_decay internally. Only a
        % tiny nonzero floor is enforced on extra_decay (numerical
        % stability against exact division-by-equal-tau, not a physical
        % constraint) -- the optimizer is free to choose anywhere from
        % a near-symmetric (single-exponential-like) shape, appropriate
        % for odors dominated by a single slow surface-adsorption/
        % clearing timescale, up to a strongly asymmetric fast-on/
        % slow-off shape, appropriate for most other odors. Forcing a
        % meaningful minimum gap (previous version) fixed odor 9's
        % shape at the cost of degrading every other odor that
        % genuinely wants a large rise/decay asymmetry.
        nParams = 5;

        lb1 = [0,     0.005, 0.001, 0.05, 0]';
        ub1 = [5,     1,     5,     20,   0.5]';  % tau_deplete large ~ negligible depletion (recovers simple monotonic rise+plateau)

        seedList = [ ...
            1, 0.05, 0.25, 0.5,  0;   % fast onset, fast depletion (strong sag)
            1, 0.05, 0.25, 2.0,  0;   % fast onset, slow depletion (mild sag)
            1, 0.05, 0.25, 20,   0;   % fast onset, ~no depletion (monotonic plateau)
            1, 0.15, 0.35, 1.0,  0;   % slower onset, moderate depletion
            1, 0.30, 0.70, 20,   0;   % very slow onset, ~no depletion
            1, 0.50, 1.00, 20,   0;   % extremely slow onset, ~no depletion
            1, 0.30, 0.70, 3,    0;   % slow onset, some depletion
            1, 0.10, 1.50, 20,   0;   % fast-ish onset, strongly slower offset (big rise/decay gap)
            1, 0.30, 0.02, 20,   0;   % slow, near-symmetric single-timescale (odor-9-like: slow surface adsorption/clearing)
            ]';

        model_fit = @pid_depletion_out;
        options = optimoptions('lsqnonlin', ...
            'MaxFunctionEvaluations',1e5, 'MaxIterations',1e5, ...
            'Display','off');

        paramsout = zeros(nParams, nOdors);
        resnormPerOdor = zeros(1, nOdors);
        exitflag = zeros(1, nOdors);
        output = cell(1, nOdors);
        residual = cell(1, nOdors);

        % Diagnostic storage: every seed's rn BEFORE optimization (at
        % the raw seed values) and AFTER optimization, for every odor.
        % This lets us tell apart two very different failure modes for
        % a given seed without needing to manually re-run anything:
        %   (a) rn_after >= rn_before  -> the optimizer failed to
        %       improve from that starting point at all -- a red flag
        %       suggesting a bug in the residual/model function, not
        %       just an unlucky local minimum.
        %   (b) rn_after < rn_before, but still worse than the winning
        %       seed's result -> a genuine local minimum: that region
        %       of parameter space really is worse for this odor, the
        %       model/optimizer are working correctly, it's just not
        %       the best available fit within this model class.
        rnBeforeAll = nan(size(seedList,2), nOdors);
        rnAfterAll  = nan(size(seedList,2), nOdors);

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
                P0 = min(max(P0, lb1), ub1);

                rnBeforeAll(s,k) = sum(depletion_residuals_single(P0, k).^2);

                try
                    [Pfit, rn, res, ef, out] = lsqnonlin( ...
                        @(P) depletion_residuals_single(P, k), P0, lb1, ub1, options);
                catch ME
                    warning('Odor %d, seed %d failed: %s', k, s, ME.message);
                    continue
                end
                rnAfterAll(s,k) = rn;

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

        % --- Automatic diagnostic printout ---
        fprintf('\n--- Depletion-model fit diagnostic (per odor, per seed) ---\n');
        for k = 1:nOdors
            if all(isnan(rnAfterAll(:,k)))
                continue
            end
            [bestRn, bestS] = min(rnAfterAll(:,k));
            fprintf('Odor %2d: winning seed = %d, rn = %8.4g\n', k, bestS, bestRn);
            for s = 1:size(seedList,2)
                if isnan(rnAfterAll(s,k))
                    continue
                end
                flag = '';
                if rnAfterAll(s,k) >= rnBeforeAll(s,k)
                    flag = '  <-- optimizer did NOT improve from this seed (possible bug)';
                elseif s ~= bestS && rnAfterAll(s,k) > 1.5*bestRn
                    flag = '  (converged to a clearly worse local minimum)';
                end
                fprintf('    seed %d: rn_before = %8.4g -> rn_after = %8.4g%s\n', ...
                    s, rnBeforeAll(s,k), rnAfterAll(s,k), flag);
            end
        end
        fprintf('--- end diagnostic ---\n\n');

        [kernelsout,PIDOut,DriveOut] = makeDepletionKernel(ValveVector,paramsout,timestep,TimeVector);

    case 'depletion1exp'
        % Same source-depletion drive as 'depletion', but with a
        % SINGLE-exponential (first-order low-pass) kernel instead of a
        % difference-of-exponentials. Physical motivation: if the slow
        % rise is really due to odor slowly accumulating in dead volume
        % (tubing etc.) that also clears slowly after valve-off, that's
        % a single time-constant process -- h(t) = (A/tau)*exp(-t/tau)
        % -- not a rise-then-decay "hump" kernel. A difference-of-
        % exponentials kernel forces tau_rise and tau_decay apart to
        % avoid degenerating toward zero, which fights against exactly
        % this symmetric-single-tau physical picture. This kernel has
        % NO such constraint: it's inherently monotonic-decaying and
        % naturally produces a slow, symmetric-ish rise-and-linger when
        % convolved with a boxcar, governed by one shared time constant.
        %
        % params per odor: [A, tau, tau_deplete, delay]
        nParams = 4;

        lb1 = [0,    0.01, 0.05, 0]';
        ub1 = [10,   2,    20,   0.5]';  % tau_deplete large ~ negligible depletion

        seedList = [ ...
            1,  0.05, 0.5,  0;   % fast clearing, fast depletion
            1,  0.05, 20,   0;   % fast clearing, ~no depletion
            1,  0.30, 20,   0;   % slow clearing (dead-volume story), ~no depletion
            1,  0.30, 3,    0;   % slow clearing, some depletion
            1,  0.60, 20,   0;   % very slow clearing, ~no depletion
            1,  0.15, 1.0,  0;   % moderate clearing, moderate depletion
            ]';

        model_fit = @pid_depletion1exp_out;
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
                P0 = min(max(P0, lb1), ub1);

                try
                    [Pfit, rn, res, ef, out] = lsqnonlin( ...
                        @(P) depletion1exp_residuals_single(P, k), P0, lb1, ub1, options);
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

        [kernelsout,PIDOut,DriveOut] = makeDepletion1expKernel(ValveVector,paramsout,timestep,TimeVector);

    case 'depletion_auto'
        % Some odors are well described by an asymmetric fast-on/
        % slow-off kernel (difference of exponentials); others --
        % dominated by a single slow surface-adsorption/clearing
        % timescale, giving a near-symmetric slow rise and slow linger
        % -- are better described by a single-exponential kernel.
        % Forcing every odor through ONE shared functional form was
        % hurting whichever group didn't match it (diff-of-exp left a
        % couple of slow-rise odors stuck in a degenerate near-zero-
        % amplitude local minimum; single-exp broke the asymmetry that
        % most other odors genuinely have). Here we fit BOTH forms per
        % odor, each with its own multi-start, and keep whichever
        % achieves the lower residual for that specific odor.

        optionsAuto = optimoptions('lsqnonlin', ...
            'MaxFunctionEvaluations',1e5, 'MaxIterations',1e5, ...
            'Display','off');

        % --- candidate A: difference-of-exponentials ---
        % params: [A, tau_rise, extra_decay, tau_deplete, delay]
        lbA = [0,    0.005, 0.001, 0.05, 0]';
        ubA = [5,    1,     5,     20,   0.5]';
        seedsA = [ ...
            1, 0.05, 0.25, 0.5,  0;
            1, 0.05, 0.25, 2.0,  0;
            1, 0.05, 0.25, 20,   0;
            1, 0.15, 0.35, 1.0,  0;
            1, 0.30, 0.70, 20,   0;
            1, 0.50, 1.00, 20,   0;
            1, 0.30, 0.70, 3,    0;
            1, 0.10, 1.50, 20,   0;
            1, 0.30, 0.02, 20,   0;
            ]';

        % --- candidate B: single-exponential ---
        % params: [A, tau, tau_deplete, delay]
        lbB = [0,    0.01, 0.05, 0]';
        ubB = [10,   2,    20,   0.5]';
        seedsB = [ ...
            1,  0.05, 0.5,  0;
            1,  0.05, 20,   0;
            1,  0.30, 20,   0;
            1,  0.30, 3,    0;
            1,  0.60, 20,   0;
            1,  0.15, 1.0,  0;
            ]';

        paramsoutA = zeros(5, nOdors);  % NaN-padded if candidate B wins
        paramsoutB = zeros(4, nOdors);
        modelUsed  = zeros(1, nOdors);  % 1 = diff-of-exp, 2 = single-exp
        resnormPerOdor = zeros(1, nOdors);
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

            bestRn = Inf; bestModel = 0; bestP = [];

            for s = 1:size(seedsA,2)
                P0 = min(max(seedsA(:,s), lbA), ubA);
                try
                    [Pfit, rn] = lsqnonlin(@(P) depletion_residuals_single(P, k), P0, lbA, ubA, optionsAuto);
                catch ME
                    warning('Odor %d, model A, seed %d failed: %s', k, s, ME.message);
                    continue
                end
                if rn < bestRn
                    bestRn = rn; bestModel = 1; bestP = Pfit;
                end
            end

            for s = 1:size(seedsB,2)
                P0 = min(max(seedsB(:,s), lbB), ubB);
                try
                    [Pfit, rn] = lsqnonlin(@(P) depletion1exp_residuals_single(P, k), P0, lbB, ubB, optionsAuto);
                catch ME
                    warning('Odor %d, model B, seed %d failed: %s', k, s, ME.message);
                    continue
                end
                if rn < bestRn
                    bestRn = rn; bestModel = 2; bestP = Pfit;
                end
            end

            resnormPerOdor(k) = bestRn;
            modelUsed(k) = bestModel;
            if bestModel == 1
                paramsoutA(:,k) = bestP;
                [K,PO,DO] = makeDepletionKernel(ValveVector, bestP, timestep, TimeVector);
            elseif bestModel == 2
                paramsoutB(:,k) = bestP;
                [K,PO,DO] = makeDepletion1expKernel(ValveVector, bestP, timestep, TimeVector);
            else
                continue % both candidates failed for this odor
            end
            kernelsout(:,k) = K;
            PIDOut(:,k)     = PO;
            DriveOut(:,k)   = DO;
        end

        resnorm  = sum(resnormPerOdor);
        exitflag = modelUsed;  % repurposed here to report which model won per odor (1=diff-exp, 2=single-exp)
        residual = PIDOut - ydata;
        output   = struct('modelUsed', modelUsed, 'paramsoutA_diffexp', paramsoutA, 'paramsoutB_singleexp', paramsoutB);
        paramsout = output; % for convenience/backward-compat when inspecting after a run

    case 'depletion_signed'
        % Generalization of 'depletion': instead of a drive that can
        % only ever DEPLETE (fall from 1 toward 0 while the valve is
        % open), let the source availability trajectory go either
        % direction. drive(t) = valve(t) * [Ainf + (1-Ainf)*exp(-t/tau_source)]
        % starts at 1 the instant the valve opens and relaxes toward
        % Ainf:
        %   Ainf < 1  -> depletion (headspace runs low during the pulse)
        %   Ainf > 1  -> "anti-depletion"/buildup (source takes time to
        %                reach full strength -- e.g. vapor slowly
        %                equilibrating in the headspace right as flow
        %                starts) -- this can produce a genuinely slow,
        %                well-behaved rise WITHOUT needing the kernel
        %                itself to degenerate toward tau_rise->0.
        %   Ainf = 1  -> no source dynamics (recovers simple monotonic
        %                rise+plateau)
        % Kernel stays the simple, always-positive, well-behaved
        % difference-of-exponentials form -- same as 'depletion'.
        %
        % params per odor: [A, tau_rise, extra_decay, Ainf, tau_source, delay]
        nParams = 6;

        lb1 = [0,    0.005, 0.001, 0,    0.05, 0]';
        ub1 = [5,    1,     5,     3,    20,   0.5]';  % Ainf in [0,3]: 0-1 = depletion, 1 = none, 1-3 = buildup

        seedList = [ ...
            1, 0.05, 0.25, 0.3,  0.5,  0;   % fast onset, depletion
            1, 0.05, 0.25, 1.0,  20,   0;   % fast onset, no source dynamics
            1, 0.05, 0.25, 2.0,  0.5,  0;   % fast onset, strong buildup (fast)
            1, 0.05, 0.25, 2.0,  2.0,  0;   % fast onset, strong buildup (slow)
            1, 0.10, 0.35, 3.0,  1.0,  0;   % fast onset, very strong buildup
            1, 0.30, 0.70, 1.0,  20,   0;   % slow onset (kernel-driven), no source dynamics
            1, 0.10, 0.35, 1.5,  0.8,  0;   % moderate buildup
            1, 0.01, 1.50, 0.3,  0.05, 0;   % fast source transient + genuinely long kernel tail (fast jump, slow sustained climb)
            1, 0.01, 2.50, 0.5,  0.05, 0;   % same, even longer tail
            ]';

        model_fit = @pid_depletion_signed_out;
        options = optimoptions('lsqnonlin', ...
            'MaxFunctionEvaluations',1e5, 'MaxIterations',1e5, ...
            'Display','off');

        paramsout = zeros(nParams, nOdors);
        resnormPerOdor = zeros(1, nOdors);
        exitflag = zeros(1, nOdors);
        output = cell(1, nOdors);
        residual = cell(1, nOdors);

        % Diagnostic: track every seed's rn before (at raw seed values)
        % and after optimization, per odor, plus which seed number won.
        rnBeforeAll = nan(size(seedList,2), nOdors);
        rnAfterAll  = nan(size(seedList,2), nOdors);
        winningSeed = zeros(1, nOdors);
        allFits     = cell(size(seedList,2), nOdors);  % store every seed's converged params, not just the winner's

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
                P0 = min(max(P0, lb1), ub1);

                rnBeforeAll(s,k) = sum(depletion_signed_residuals_single(P0, k).^2);

                try
                    [Pfit, rn, res, ef, out] = lsqnonlin( ...
                        @(P) depletion_signed_residuals_single(P, k), P0, lb1, ub1, options);
                catch ME
                    warning('Odor %d, seed %d failed: %s', k, s, ME.message);
                    continue
                end
                rnAfterAll(s,k) = rn;
                allFits{s,k}    = Pfit;

                if rn < bestResnorm
                    bestResnorm       = rn;
                    paramsout(:,k)    = Pfit;
                    resnormPerOdor(k) = rn;
                    exitflag(k)       = ef;
                    output{k}         = out;
                    residual{k}       = res;
                    winningSeed(k)    = s;
                end
            end
        end
        resnorm = sum(resnormPerOdor);

        fprintf('\n--- depletion_signed: per-odor summary of new seeds (8,9) vs winner ---\n');
        for k = 1:nOdors
            if isnan(rnAfterAll(1,k)), continue; end
            bestOfOriginal7 = min(rnAfterAll(1:7,k));
            bestOfNew2      = min(rnAfterAll(8:9,k));
            newSeedWonStr = '';
            if bestOfNew2 < bestOfOriginal7
                newSeedWonStr = '  <-- new seed WON here';
            end
            fprintf('Odor %2d: winning seed = %d (rn=%.6g) | best of original 7 = %.6g | best of new seeds (8,9) = %.6g%s\n', ...
                k, winningSeed(k), resnormPerOdor(k), bestOfOriginal7, bestOfNew2, newSeedWonStr);
        end
        fprintf('--- end summary ---\n\n');

        fprintf('\n--- depletion_signed: per-seed diagnostic (odor 9 shown; edit odorsToShow to see others) ---\n');
        odorsToShow = 9;
        for k = odorsToShow
            fprintf('Odor %2d: winning seed = %d, rn = %.10g\n', k, winningSeed(k), resnormPerOdor(k));
            for s = 1:size(seedList,2)
                if isnan(rnAfterAll(s,k)), continue; end
                marker = '';
                if s == winningSeed(k), marker = '  <-- WINNER'; end
                fprintf('    seed %2d: rn_before = %12.6g -> rn_after = %.10g%s\n', ...
                    s, rnBeforeAll(s,k), rnAfterAll(s,k), marker);
                fprintf('              params = [%s]\n', sprintf('%.6g  ', allFits{s,k}));
            end
        end
        fprintf('--- end diagnostic ---\n\n');

        fprintf('\n--- depletion_signed: Ainf per odor (< 1 = depletion, > 1 = buildup, = 1 = neither) ---\n');
        for k = 1:nOdors
            if all(paramsout(:,k)==0), continue; end
            fprintf('Odor %2d: Ainf = %.3f, tau_source = %.3f, tau_rise = %.4f, tau_decay = %.4f\n', ...
                k, paramsout(4,k), paramsout(5,k), paramsout(2,k), paramsout(2,k)+paramsout(3,k));
        end
        fprintf('--- end summary ---\n\n');

        [kernelsout,PIDOut,DriveOut] = makeDepletionSignedKernel(ValveVector,paramsout,timestep,TimeVector);

    case 'stretched_test'
        % Full run: STRETCHED-exponential decay tail
        % (exp(-(t/tau_decay)^beta), beta<1 gives a fatter/slower tail
        % than a plain exponential) across all odors. beta=1 recovers
        % the plain-exponential 'depletion_signed' model exactly, so
        % this is a strict generalization -- for each odor we can see
        % directly whether it wants beta<1 (genuinely non-exponential
        % relaxation) or converges back to beta~1 (plain exponential
        % was already fine for that odor).
        %
        % params: [A, tau_rise, tau_decay, beta, Ainf, tau_source, delay]
        targetOdors = 1:nOdors;

        nParams = 7;
        lb1 = [0,    0.005, 0.02, 0.2, 0,   0.05, 0]';
        ub1 = [5,    1,     3,    1.0, 3,   20,   0.5]';  % beta=1 -> plain exponential

        seedList = [ ...
            1,    0.05, 0.30, 1.0, 1.0, 20,   0;      % generic: fast onset, no source dynamics, plain exponential
            1,    0.05, 0.30, 0.5, 0.3, 0.5,  0;       % generic: fast onset, depletion, moderate stretch
            1,    0.05, 0.30, 0.5, 1.0, 20,   0;       % generic: fast onset, no source dynamics, stretched
            1,    0.10, 0.50, 0.5, 2.0, 1.0,  0;       % generic: buildup, stretched
            1,    0.30, 0.80, 0.6, 1.0, 20,   0;       % generic: slow onset (kernel-driven), no source dynamics
            0.10, 0.01, 0.17, 1.0, 0.57, 0.05, 0.16;   % odor-9-like: fast jump + fast depletion, plain exponential
            0.10, 0.01, 0.30, 0.5, 0.57, 0.05, 0.16;   % odor-9-like: fast jump + fast depletion, stretched
            0.10, 0.01, 0.80, 0.35,0.57, 0.05, 0.16;   % odor-9-like: fast jump + fast depletion, strongly stretched
            ]';

        options = optimoptions('lsqnonlin', ...
            'MaxFunctionEvaluations',1e5, 'MaxIterations',1e5, ...
            'Display','off');

        paramsout = zeros(nParams, nOdors);
        resnormPerOdor = nan(1, nOdors);
        exitflag = zeros(1, nOdors);
        output = cell(1, nOdors);
        residual = cell(1, nOdors);
        kernelsout = zeros(size(xdata,1), nOdors);
        PIDOut     = zeros(size(xdata,1), nOdors);
        DriveOut   = zeros(size(xdata,1), nOdors);

        % previous plain-exponential ('depletion_signed') resnorm per
        % odor, for direct before/after comparison -- from the
        % confirmed run earlier in this conversation.
        prevRn = [277.759, 1.11356, 0.184535, 0.616131, 0.294279, 0.253277, ...
                  0.340888, 0.905693, 0.828077, 1266.29, 0.181477, 0.341709, ...
                  1.05098, 0.306224, 1.36891, 22.1813];

        fprintf('\n--- stretched_test: all odors, beta<1 vs plain-exponential baseline ---\n');
        for k = targetOdors
            bestResnorm = Inf; bestBeta = NaN;
            for s = 1:size(seedList,2)
                P0 = min(max(seedList(:,s), lb1), ub1);
                try
                    [Pfit, rn] = lsqnonlin(@(P) stretched_residuals_single(P, k), P0, lb1, ub1, options);
                catch ME
                    warning('Odor %d, seed %d failed: %s', k, s, ME.message);
                    continue
                end
                if rn < bestResnorm
                    bestResnorm = rn; paramsout(:,k) = Pfit; bestBeta = Pfit(4);
                end
            end
            resnormPerOdor(k) = bestResnorm;

            pctImprove = 100*(prevRn(k) - bestResnorm)/prevRn(k);
            verdict = '';
            if pctImprove > 5
                verdict = '  <-- MEANINGFUL improvement';
            elseif pctImprove < -1
                verdict = '  <-- WORSE (unexpected, check seeds)';
            end
            fprintf('Odor %2d: beta=%.3f | rn: %.5g -> %.5g (%.1f%%)%s\n', ...
                k, bestBeta, prevRn(k), bestResnorm, pctImprove, verdict);

            [K,PO,DO] = makeStretchedKernel(ValveVector, paramsout(:,k), timestep, TimeVector);
            kernelsout(:,k) = K; PIDOut(:,k) = PO; DriveOut(:,k) = DO;
        end
        fprintf('--- end stretched_test ---\n\n');

        resnorm = sum(resnormPerOdor(~isnan(resnormPerOdor)));

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

function [res] = depletion_residuals_single(P, k)
    % Residual vector for lsqnonlin, single odor k.
    zdata = pid_depletion_out(P, xdata(:,k));
    res = zdata - ydata(:,k);
end

function [zdata] = pid_depletion_out(P, xdata_col)
    % Source-depletion model: drive signal decays exponentially while
    % the valve is open (headspace depletion), convolved with a simple
    % always-positive kernel (fast rise, slower decay).
    A           = P(1);
    tau_rise    = P(2);
    extra_decay = P(3);
    tau_deplete = P(4);
    delay       = P(5);
    tau_decay   = tau_rise + extra_decay; % ensures tau_decay > tau_rise by a meaningful margin

    N  = size(xdata_col,1);
    dt = timestep;                  % captured from outer scope
    tvec = (0:N-1)' * dt;

    % TimeVector (captured from outer scope) already references t=0 at
    % odor onset, so it doubles directly as "time since onset" for the
    % depletion clock. Multiplying by xdata_col (the 0/1 valve trace)
    % ensures depletion only accumulates while the valve is open, and
    % is exactly 0 when the valve is closed (matching the raw drive).
    ageSinceOnset   = TimeVector(:);
    depletionFactor = exp(-max(ageSinceOnset,0)/tau_deplete);
    effectiveDrive  = xdata_col .* depletionFactor;

    t = tvec - delay;
    h = A * (exp(-t/tau_decay) - exp(-t/tau_rise));
    h(t < 0) = 0;   % causal kernel

    full = conv(effectiveDrive', h', 'full');
    zdata = full(1:N)';
end

function [K,PIDOut,DriveOut] = makeDepletionKernel(ValveVector,Params,timestep,TimeVectorIn)
    N  = numel(ValveVector);
    dt = timestep;
    tvec = (0:N-1)' * dt;
    nOd = size(Params,2);

    K        = zeros(N,nOd);
    PIDOut   = zeros(N,nOd);
    DriveOut = zeros(N,nOd);
    ageSinceOnset = TimeVectorIn(:);

    for kk = 1:nOd
        A           = Params(1,kk);
        tau_rise    = Params(2,kk);
        extra_decay = Params(3,kk);
        tau_deplete = Params(4,kk);
        delay       = Params(5,kk);
        tau_decay   = tau_rise + extra_decay; % ensures tau_decay > tau_rise by a meaningful margin

        depletionFactor = exp(-max(ageSinceOnset,0)/tau_deplete);
        drive = ValveVector(:) .* depletionFactor;

        t = tvec - delay;
        h = A * (exp(-t/tau_decay) - exp(-t/tau_rise));
        h(t < 0) = 0;

        full = conv(drive', h', 'full');
        PIDOut(:,kk)   = full(1:N)';
        K(:,kk)        = h;
        DriveOut(:,kk) = drive;
    end
end

function [res] = depletion1exp_residuals_single(P, k)
    zdata = pid_depletion1exp_out(P, xdata(:,k));
    res = zdata - ydata(:,k);
end

function [zdata] = pid_depletion1exp_out(P, xdata_col)
    % Source-depletion drive (same as pid_depletion_out), convolved with
    % a SINGLE-exponential (first-order low-pass) kernel: h(t) =
    % (A/tau)*exp(-(t-delay)/tau), monotonic decay, no rise/decay
    % time-constant split -- one shared timescale governs both the
    % accumulation (rise) and the clearing (offset linger).
    A           = P(1);
    tau         = P(2);
    tau_deplete = P(3);
    delay       = P(4);

    N  = size(xdata_col,1);
    dt = timestep;              % captured from outer scope
    tvec = (0:N-1)' * dt;

    ageSinceOnset   = TimeVector(:);
    depletionFactor = exp(-max(ageSinceOnset,0)/tau_deplete);
    effectiveDrive  = xdata_col .* depletionFactor;

    t = tvec - delay;
    h = (A/tau) * exp(-t/tau);
    h(t < 0) = 0;   % causal kernel

    full = conv(effectiveDrive', h', 'full');
    zdata = full(1:N)';
end

function [K,PIDOut,DriveOut] = makeDepletion1expKernel(ValveVector,Params,timestep,TimeVectorIn)
    N  = numel(ValveVector);
    dt = timestep;
    tvec = (0:N-1)' * dt;
    nOd = size(Params,2);

    K        = zeros(N,nOd);
    PIDOut   = zeros(N,nOd);
    DriveOut = zeros(N,nOd);
    ageSinceOnset = TimeVectorIn(:);

    for kk = 1:nOd
        A           = Params(1,kk);
        tau         = Params(2,kk);
        tau_deplete = Params(3,kk);
        delay       = Params(4,kk);

        depletionFactor = exp(-max(ageSinceOnset,0)/tau_deplete);
        drive = ValveVector(:) .* depletionFactor;

        t = tvec - delay;
        h = (A/tau) * exp(-t/tau);
        h(t < 0) = 0;

        full = conv(drive', h', 'full');
        PIDOut(:,kk)   = full(1:N)';
        K(:,kk)        = h;
        DriveOut(:,kk) = drive;
    end
end

function [res] = depletion_signed_residuals_single(P, k)
    zdata = pid_depletion_signed_out(P, xdata(:,k));
    res = zdata - ydata(:,k);
end

function [zdata] = pid_depletion_signed_out(P, xdata_col)
    % Generalized source-dynamics drive (can deplete OR build up),
    % convolved with the same well-behaved, always-positive
    % difference-of-exponentials kernel used in the plain 'depletion'
    % case.
    A           = P(1);
    tau_rise    = P(2);
    extra_decay = P(3);
    Ainf        = P(4);
    tau_source  = P(5);
    delay       = P(6);
    tau_decay   = tau_rise + extra_decay;

    N  = size(xdata_col,1);
    dt = timestep;              % captured from outer scope
    tvec = (0:N-1)' * dt;

    ageSinceOnset = TimeVector(:);
    sourceFactor  = Ainf + (1-Ainf) * exp(-max(ageSinceOnset,0)/tau_source);
    effectiveDrive = xdata_col .* sourceFactor;

    t = tvec - delay;
    h = A * (exp(-t/tau_decay) - exp(-t/tau_rise));
    h(t < 0) = 0;   % causal kernel

    full = conv(effectiveDrive', h', 'full');
    zdata = full(1:N)';
end

function [K,PIDOut,DriveOut] = makeDepletionSignedKernel(ValveVector,Params,timestep,TimeVectorIn)
    N  = numel(ValveVector);
    dt = timestep;
    tvec = (0:N-1)' * dt;
    nOd = size(Params,2);

    K        = zeros(N,nOd);
    PIDOut   = zeros(N,nOd);
    DriveOut = zeros(N,nOd);
    ageSinceOnset = TimeVectorIn(:);

    for kk = 1:nOd
        A           = Params(1,kk);
        tau_rise    = Params(2,kk);
        extra_decay = Params(3,kk);
        Ainf        = Params(4,kk);
        tau_source  = Params(5,kk);
        delay       = Params(6,kk);
        tau_decay   = tau_rise + extra_decay;

        sourceFactor = Ainf + (1-Ainf) * exp(-max(ageSinceOnset,0)/tau_source);
        drive = ValveVector(:) .* sourceFactor;

        t = tvec - delay;
        h = A * (exp(-t/tau_decay) - exp(-t/tau_rise));
        h(t < 0) = 0;

        full = conv(drive', h', 'full');
        PIDOut(:,kk)   = full(1:N)';
        K(:,kk)        = h;
        DriveOut(:,kk) = drive;
    end
end

function [res] = stretched_residuals_single(P, k)
    zdata = pid_stretched_out(P, xdata(:,k));
    res = zdata - ydata(:,k);
end

function [zdata] = pid_stretched_out(P, xdata_col)
    % Same signed source-dynamics drive as pid_depletion_signed_out, but
    % with a STRETCHED-exponential decay term in the kernel:
    % exp(-(t/tau_decay)^beta). beta=1 recovers a plain exponential
    % (identical to pid_depletion_signed_out); beta<1 gives a fatter,
    % slower-converging tail -- less front-loaded curvature approaching
    % the plateau.
    A          = P(1);
    tau_rise   = P(2);
    tau_decay  = P(3);
    beta       = P(4);
    Ainf       = P(5);
    tau_source = P(6);
    delay      = P(7);

    N  = size(xdata_col,1);
    dt = timestep;              % captured from outer scope
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

function [K,PIDOut,DriveOut] = makeStretchedKernel(ValveVector,P,timestep,TimeVectorIn)
    N  = numel(ValveVector);
    dt = timestep;
    tvec = (0:N-1)' * dt;
    ageSinceOnset = TimeVectorIn(:);

    A          = P(1);
    tau_rise   = P(2);
    tau_decay  = P(3);
    beta       = P(4);
    Ainf       = P(5);
    tau_source = P(6);
    delay      = P(7);

    sourceFactor = Ainf + (1-Ainf) * exp(-max(ageSinceOnset,0)/tau_source);
    DriveOut = ValveVector(:) .* sourceFactor;

    t = tvec - delay;
    t(t < 0) = NaN;
    K = A * (exp(-(t/tau_decay).^beta) - exp(-t/tau_rise));
    K(isnan(K)) = 0;

    full = conv(DriveOut', K', 'full');
    PIDOut = full(1:N)';
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