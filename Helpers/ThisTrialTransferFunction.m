function [TF] = ThisTrialTransferFunction(target)

TF_bins = 100;

lever_max = 4.8; %triggerUpLim
lever_min = 0.2; %triggerLowLim
zone_width = 0.3;
gain = 1;
total_motor_locations = 100;
minimumtarget = 1;

% calculate stepsize - lever displacement corresponding to one location
stepsize = (lever_max - minimumtarget)/(total_motor_locations + 0.5);
% compute number of locations to be allocated to the target zone
locations_per_zone = round(zone_width/stepsize);

% rescale stepsize if needed for gain perturbation trials
stepsize = stepsize*gain;

start_location = numel(target:stepsize:lever_max);
end_location = -numel(target:-stepsize:lever_min);
TF = linspace(end_location,start_location,TF_bins);

% safetychecks
TF = round(TF);
TF(TF>115) = 115;
TF(TF<-115) = -115;



end