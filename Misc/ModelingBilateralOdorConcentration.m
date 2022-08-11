Locations = [-100:1:100];
mu = 0; sigma = 30;

% each unit in odor space = 0.5 mm; 100 units = 50 mm
% inter-nostril distance = 2 mm, 4 units in odor space
delta_nostrils = 2;

C_center = normpdf(Locations,mu,sigma);
C_center = C_center/max(C_center);
C_left = normpdf(Locations+delta_nostrils,mu,sigma);
C_left = C_left/max(C_left);
C_right = normpdf(Locations-delta_nostrils,mu,sigma);
C_right = C_right/max(C_right);

figure, 
subplot(1,2,1);
plot(Locations, C_center, 'k');
hold on
plot(Locations, C_left, 'b');
plot(Locations, C_right, 'r');

subplot(1,2,2);
plot(Locations, C_left - C_right, 'k');


