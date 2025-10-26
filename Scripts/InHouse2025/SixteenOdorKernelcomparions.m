load('/mnt/data/sniffGLM/T2_20250521_r0_shared/SixteenOdorsFitOutput.mat');
figure;

AirKernels = squeeze(kernels.pGLM(:,1,:))';
OdorKernels = kernels.pGLM(:,2:end,:);
OdorKernels = reshape(OdorKernels,70,[])';

subplot(1,6,1);
imagesc(AirKernels,[-5 5]);
%set(gca, 'YLim', [0 150])

subplot(1,6,3);
imagesc(OdorKernels,[-5 5]);
%set(gca, 'YLim', [0 2400])

subplot(1,6,5);
f = [];
for i = 1:size(OdorKernels,1)
    [~,f(i)] = max(abs(OdorKernels(i,:)));
end
[~,S] = sort(f);
imagesc(OdorKernels(S,:), [-5 5]);
%set(gca, 'YLim', [0 2400])

load('/mnt/data/sniffGLM/T3_20250516_r0_shared/SixteenOdorsFitOutput.mat');

AirKernels = squeeze(kernels.pGLM(:,1,:))';
OdorKernels = kernels.pGLM(:,2:end,:);
OdorKernels = reshape(OdorKernels,70,[])';

subplot(1,6,2);
imagesc(AirKernels, [-5 5]);
%set(gca, 'YLim', [0 150])

subplot(1,6,4);
imagesc(OdorKernels, [-5 5]);
%set(gca, 'YLim', [0 2400])

subplot(1,6,6);
f = [];
for i = 1:size(OdorKernels,1)
    [~,f(i)] = max(abs(OdorKernels(i,:)));
end
[~,S] = sort(f);
imagesc(OdorKernels(S,:), [-5 5]);
%set(gca, 'YLim', [0 2400])


colormap(brewermap([100],'RdBu'))




