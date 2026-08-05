figure;
for i = 1:16
    subplot(4,4,i);
    kernelparams = F.paramsout(:,i);
    A            = kernelparams(1); % ON kernel scakar
tau_rise     = kernelparams(2); % ON kernel fast tau
tau_decay_on = kernelparams(3); % ON kernel slow tau
beta_on      = kernelparams(4); % ON kernel slow modifier

    dt = 0.01; % 10 ms
    N = 0.5/dt; % 1 sec
    tvec = (0:N-1)' * dt;
    fastOne = exp(-tvec/tau_rise);
    slowOne = exp(-(tvec/tau_decay_on).^beta_on);
    slowTwo = exp(-(tvec/tau_decay_on));
    plot(tvec,fastOne,'r');
    hold on
    plot(tvec,slowOne,'k');
    plot(tvec,slowTwo,'b');

    title([num2str(tau_rise),' ',num2str(tau_decay_on),' ',num2str(beta_on)])
end

odorNames = {'PAnisD', 'Hept', 'ET', 'AlTig', ...
    'EtProp', 'EtBut', 'EtVal', 'IAA', ...
    'PropTig', 'PhEtAc', 'PropBut', 'Hex', ...
    'AcPh', 'ValDe','GamTurp', 'Oil'};

figure; 
for p = 1:10
    subplot(2,5,p);
    plot(1:16, F.paramsout(p,:), '.', 'MarkerSize',10);
    set(gca,'XTick',[1:16],'XTickLabel',odorNames, 'XTickLabelRotation', 45, 'XGrid', 'on');
end
