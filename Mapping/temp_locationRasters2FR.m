[PassiveTuning] = LocationTuningDataParser('E3_20220713_o3.mat');

% some info for the raster size
% window = [-TuningParams(1,8)/2 (sum(TuningParams(4:7)) + TuningParams(1,8)/2)]/1000;
% PassiveTuning.FileLocations = FileLocations;
% PassiveTuning.ITI = 1000*[-window(1) TuningParams(1,8)/1000+window(1)];
% PassiveTuning.PreOdor = TuningParams(1,4);
% PassiveTuning.Odor = TuningParams(1,5);
% PassiveTuning.PostOdor = TuningParams(1,6) + TuningParams(1,7);

% chop off the rasters to have shorter pre-odor and post-odor stretches
Air1 = (PassiveTuning.ITI(1) + PassiveTuning.PreOdor) + [(1-PassiveTuning.Odor) 0];
Odor = (PassiveTuning.ITI(1) + PassiveTuning.PreOdor) + [1 PassiveTuning.Odor];
Air2 = Odor(2) + [1 + (1+diff(Odor))];

% convert Spike rasters to FR
timewindow = 200; % time window over which smoothing occurs - change as needed
t_wid = timewindow;  % width of kernel (in ms)
taxis = -(t_wid*5):(t_wid*5);  % e.g. make a time axis of 1000 ms for a window of 100 ms
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

%% PSTH
[Nneurons, Nodor, Nloc, lastTS, Nrep] = size(PassiveTuning.RasterOut);

%% Smooth PSTH - to get a continuous (time-dependent) rate variable

clusterNum = 1:Nneurons;
for clusterIdx = 1:size(PassiveTuning.RasterOut,1) % each unit
    for x = 1:size(PassiveTuning.RasterOut,2) % each odor
        for y = 1:size(PassiveTuning.RasterOut,3) % each location
            for j = 1:size(PassiveTuning.RasterOut,5) % each repeat
                tempPSTH = squeeze(PassiveTuning.RasterOut(clusterIdx,x,y,:,j)); % for one cluster, all times for one type
                tempPSTH(:,(Air2(2)+1):end) = [];
                tempPSTH(:,1:(Air1(1)-1)) = [];
                zs = conv(tempPSTH,gauss_kernel,'same'); %in ms
                PSTH(clusterIdx,x,y,:,j) = zs*1000; %converting firing rate to Hz (1 ms = 1000 Hz)
            end
        end
    end
end

