function residuals = PSTHresiduals(StartingKernels,InputVector,PSTH)
[baseline,kernels,locationcoef] = ParseSniffKernels(StartingKernels);
[zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,InputVector);
residuals = mean((PSTH - zdata).^2);
end