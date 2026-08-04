function PIDrecon = reconstructTwoRegimeFit(DriveOut, kernelOn, offKernel, boundaryValue, switchIdx)
% Reconstructs a fitted PID trace from the pieces stored in Fitted,
% without needing to rerun any fitting code. Exactly reproduces
% Fitted.PIDOut(:,k) for a given odor k when called as:
%
%   PIDrecon = reconstructTwoRegimeFit( ...
%       Fitted.DriveOut(:,k), Fitted.kernelsout(:,k), ...
%       Fitted.offKernel(:,k), Fitted.boundaryValue(k), Fitted.switchIdx(k));
%
% This is NOT a single convolution end-to-end -- the ON regime is a
% convolution of the drive against the ON kernel, but the OFF regime is
% a separate process (a scaled copy of the OFF kernel, anchored to
% whatever the ON regime's own result was at the moment the valve
% closure is felt). The two pieces are stitched together at switchIdx,
% not summed as two parallel convolutions -- see the model's docstring
% in FitCIDPIDKernels.m for why (the OFF regime's amplitude is *derived
% from* the ON regime's result, not an independent input signal).
%
% Inputs (all column vectors/scalars for one odor):
%   DriveOut      - the depleting drive signal, Fitted.DriveOut(:,k)
%   kernelOn      - the ON-regime kernel shape, Fitted.kernelsout(:,k)
%   offKernel     - the OFF-regime kernel shape (tau=0 at index 1),
%                   Fitted.offKernel(:,k)
%   boundaryValue - scalar, the ON regime's value at the switch point,
%                   Fitted.boundaryValue(k)
%   switchIdx     - scalar sample index where the OFF regime begins,
%                   Fitted.switchIdx(k)
%
% Output:
%   PIDrecon - reconstructed trace, same length as DriveOut/kernelOn

DriveOut = DriveOut(:);
kernelOn = kernelOn(:);
offKernel = offKernel(:);
N = numel(DriveOut);

% ON regime: standard convolution
full = conv(DriveOut', kernelOn', 'full');
y_on = full(1:N)';

PIDrecon = y_on;

% OFF regime: scaled copy of the OFF kernel shape, spliced in after
% switchIdx. This is equivalent to convolving an impulse of height
% boundaryValue (at the switch sample) against offKernel -- an impulse
% convolved with any kernel is just that kernel, scaled and shifted, so
% direct placement is mathematically identical to that convolution.
if switchIdx < N
    nTail = N - switchIdx;
    PIDrecon(switchIdx+1:end) = boundaryValue * offKernel(1:nTail);
end

end