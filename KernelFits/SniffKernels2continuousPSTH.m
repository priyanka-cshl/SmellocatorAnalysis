function [zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,regressors)

if numel(locationcoef)~=3 % not independent coeffs for different odors
    locationcoef = repmat(locationcoef,1,3);
end

% regressors is 6 x n matrix, where n is no. of timebins
% rows are all sniffs, sniffs with manifold air ON, odor1 sniffs, odor2 sniffs, odor3 sniffs and Motor location 

%% model
% scale odor sniffs by the location coefficients 
for i = 1:3
    regressors(i+2,:) = regressors(i+2,:).* ...
        exp( abs(regressors(6,:)) * locationcoef(i) ) ;
end

convmode = 'full'; %'same';
zdata   = baseline + ...
    conv(regressors(1,:),kernels(:,1)',convmode) + ...
    conv(regressors(2,:),kernels(:,2)',convmode) + ...
    conv(regressors(3,:),kernels(:,3)',convmode) + ...
    conv(regressors(4,:),kernels(:,4)',convmode) + ...
    conv(regressors(5,:),kernels(:,5)',convmode) ;

if strcmp(convmode, 'full')
    bins2del = size(kernels,1) - 1;
    zdata(:,end-bins2del+1:end) = [];
%     
%     
%     % delete a few bins from the traces
%     bins2del = ceil(size(kernels,1)/2);
%     zdata(:,1:bins2del) = [];
%     zdata(:,end-bins2del+2:end) = [];
end

% ignore datapoints during perturbation
perturbationbins = find(regressors(7,:));
zdata(:,perturbationbins) = 0;

end