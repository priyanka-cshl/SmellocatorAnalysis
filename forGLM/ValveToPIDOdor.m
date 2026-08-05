function [OdorOut] = ValveToPIDOdor(valveState,kernelparams,fitdt)
% a function to get the resultant PID waveform out 
% given an input valve trace
% Inputs: valveState [t x 2]: time by Valvestate
%       : kernelparams [10 x 1];

if nargin<3
    fitdt = 0.01; % was fit at 10 ms resolution
end

%% Extract parameters
A            = kernelparams(1); % ON kernel scakar
tau_rise     = kernelparams(2); % ON kernel fast tau
tau_decay_on = kernelparams(3); % ON kernel slow tau
beta_on      = kernelparams(4); % ON kernel slow modifier

Ainf         = kernelparams(5); % Min headspace floor
tau_source   = kernelparams(6); % headspace tau

delay        = kernelparams(7); % ON dead delay

tau_decay_off_slow = kernelparams(8); % decay long tail time constant
w_slow       = kernelparams(9);  % slow decaying fraction after valve shutoff
delay_off    = kernelparams(10); % OFF dead delay

%% Make time axis and required valve traces for convolution

% available odor
valveON = find(valveState(:,2)==1,1,'first');
valveON = valveState(valveON,1); % time-point when odor valve turned on
currentT  = valveState(:,1) - valveON;
sourceFactor   = Ainf + (1-Ainf) * exp(-max(currentT,0)/tau_source);
DriveOut_ = valveState(:,2) .* sourceFactor;

% ON kernel
N  = size(valveState,1);
dt = mode(diff(valveState(:,1)));
tvec = (0:N-1)' * dt;
t = tvec - delay;
t(t < 0) = NaN;
K = A * (exp(-(t/tau_decay_on).^beta_on) - exp(-t/tau_rise));
K(isnan(K)) = 0;
K(K<0) = 0;

full = conv(DriveOut_', K', 'full');
OdorTemp = full(1:N)';

% OFF dynamics
offKernel = (1-w_slow)*exp(-tvec/tau_rise) + w_slow*exp(-tvec/tau_decay_off_slow);

% find when to switch to a non-convolutional decay dynamics
rawSwitchIdx = find(valveState(:,2)==1, 1, 'last');
OdorOut = OdorTemp;
boundaryValue = 0;
switchIdxOut  = N;
if ~isempty(rawSwitchIdx)
    switchIdxOut = min(rawSwitchIdx + round(delay_off/dt), N);
    if switchIdxOut < N
        boundaryValue = OdorTemp(switchIdxOut);
        nTail = N - switchIdxOut;
        OdorOut(switchIdxOut+1:end) = boundaryValue * offKernel(1:nTail);
    end
end

% scale for dt 
OdorOut = OdorOut/(fitdt/dt); % params were calculated at 10 ms time bin

end

