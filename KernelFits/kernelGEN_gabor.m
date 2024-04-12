function [fit_kernel] = kernelGEN_gabor(startparams,z)
t = (1:z)';
a = startparams(1);
w = startparams(2)*startparams(2);
f = startparams(3);
for x = 1:size(t)
    T = t(x,1) - startparams(5);
    Tsq = -1*(T*T);
    F = 2*pi*((f*T)+startparams(4));
    G1(t(x,1),1) = a*exp(Tsq/w)*cos(F);
end
fit_kernel = G1; % G2 is the inhibitory gaussian
end