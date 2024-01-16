% simple script to understand Modulation index from Audette et al. 2022

% settings 
N       = 10000;    % no. of random number pairs to generate 
theta   = -45;      % angle w.r.t. origin at which modulation is zero - unity line
plot_   = 1;
dot_sz  = 10;

if plot_
    figure; 
end

% generate a two column list of random numbers in the range 0 to 1
XY_pairs    = rand(N,2);
% rescale the list to go from -1 to 1
XY_pairs    = (XY_pairs - 0.5)/0.5;
if plot_
    subplot(2,3,1);
    axis square;
    scatter(XY_pairs(:,1),XY_pairs(:,2),dot_sz,[0.7 0.7 0.7],'filled');
    % highlight a few entries
    chosen_points = floor(rand(10,1)*N);
    hold on
    scatter(XY_pairs(chosen_points,1),XY_pairs(chosen_points,2),dot_sz,'k','filled');
    set(gca,'XLim',[-1.5 1.5], 'YLim',[-1.5 1.5]);
    set(gca,'XTick', [], 'YTick', []);
end

% rotate each data point to recenter on the unity line
R           = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; 
XY_rot      = (R*XY_pairs.').';
if plot_
    subplot(2,3,4);
    axis square;
    scatter(XY_rot(:,1),XY_rot(:,2),dot_sz,[0.7 0.7 0.7],'filled');
    % highlight the same entries
    % chosen_points = floor(rand(10,1)*N);
    hold on
    scatter(XY_rot(chosen_points,1),XY_rot(chosen_points,2),dot_sz,'k','filled');
    set(gca,'XLim',[-1.5 1.5], 'YLim',[-1.5 1.5]);
    set(gca,'XTick', [], 'YTick', []);
end

% calucalate the ratio (slope) for each XY pair in the rotated space
Slopes      = atand(XY_pairs(:,2)./XY_pairs(:,1));
Slopes_rot  = atand(XY_rot(:,2)./XY_rot(:,1));

if plot_
    subplot(2,3,2);
    axis square;
    colormap(brewermap([10],'RdBu'));
    scatter(XY_pairs(:,1),XY_pairs(:,2),dot_sz,Slopes,'filled');
    hold on
    scatter(XY_pairs(chosen_points,1),XY_pairs(chosen_points,2),dot_sz+5,'k','filled');
    colorbar;
    set(gca,'XTick', [], 'YTick', []);
    
    subplot(2,3,5);
    axis square;
    colormap(brewermap([10],'RdBu'));
    scatter(XY_pairs(:,1),XY_pairs(:,2),dot_sz,Slopes_rot,'filled');
    hold on
    scatter(XY_pairs(chosen_points,1),XY_pairs(chosen_points,2),dot_sz+5,'k','filled');
    colorbar;
    set(gca,'XTick', [], 'YTick', []);
end

% need to normalize the slopes by the sign of the 'X' data point to avoid
% discontinuities
Slopes      = Slopes.* (abs(XY_pairs(:,1))./XY_pairs(:,1)) / 45;
Slopes_rot  = Slopes_rot.* (abs(XY_rot(:,1))./XY_rot(:,1)) / 45;
if plot_
    subplot(2,3,3);
    axis square;
    colormap(brewermap([10],'RdBu'));
    scatter(XY_pairs(:,1),XY_pairs(:,2),dot_sz,Slopes,'filled');
    hold on
    scatter(XY_pairs(chosen_points,1),XY_pairs(chosen_points,2),dot_sz+5,'k','filled');
    colorbar;
    set(gca,'XTick', [], 'YTick', []);
    
    subplot(2,3,6);
    axis square;
    colormap(brewermap([10],'RdBu'));
    scatter(XY_pairs(:,1),XY_pairs(:,2),dot_sz,Slopes_rot,'filled');
    hold on
    scatter(XY_pairs(chosen_points,1),XY_pairs(chosen_points,2),dot_sz+5,'k','filled');
    colorbar;
    set(gca,'XTick', [], 'YTick', []);
end

