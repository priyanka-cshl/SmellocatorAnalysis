function [zdata] = SniffKernels2PSTH(baseline,kernels,locationcoef,sniffdata)

if numel(locationcoef)~=3 % not independent coeffs for different odors
    locationcoef = repmat(locationcoef,1,3);
end

x_data = sniffdata;
snifflengths    = x_data(:,11:12);

% initialize PSTH matrix
z_data = zeros(size(x_data,1),max(x_data(:,11)));

for i = 1:size(x_data,1) % every sniff
    % this sniff
    stimstate1 = x_data(i,1); % -1 for ITI, 0, 1, 2, 3 for air, odor 1-3
    location1  = exp(abs(x_data(i,2)) * locationcoef);
    
    R1 = kernels(:,1) + ...
        (stimstate1>=0)*kernels(:,2) + ...
        (stimstate1==1)*location1(1)*kernels(:,3) + ...
        (stimstate1==2)*location1(2)*kernels(:,4) + ...
        (stimstate1==3)*location1(3)*kernels(:,5);
    
    % previous sniff
    stimstate2 = x_data(i,7); % -1 for ITI, 0, 1, 2, 3 for air, odor 1-3
    location2  = exp(abs(x_data(i,8)) * locationcoef);
    
    R2 = kernels(:,1) + ...
        (stimstate2>=0)*kernels(:,2) + ...
        (stimstate2==1)*location2(1)*kernels(:,3) + ...
        (stimstate2==2)*location2(2)*kernels(:,4) + ...
        (stimstate2==3)*location2(3)*kernels(:,5);
    
    bins2delete = snifflengths(i,2);
    if length(R2)>bins2delete
        R2(1:bins2delete,:) = [];
        R2(end+1:length(R1),:) = 0;
    else
        R2 = R2*0;
    end
    
    R1 = R1 + R2 + baseline;
    z_data(i,1:length(R1)) = R1;
end

zdata = z_data;
zdata(zdata<0) = 0;
for k = 1:size(zdata,1)
    zdata(k,(1+snifflengths(k,1)):end) = 0;
end

end