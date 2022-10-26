% get spike rasters from the raw data
[PassiveTuning] = LocationTuningDataParser('E3_20220713_o3.mat');

% chop off the rasters to have shorter pre-odor and post-odor stretches
Air1 = (PassiveTuning.ITI(1) + PassiveTuning.PreOdor) + [(1-PassiveTuning.Odor) 0];
Odor = (PassiveTuning.ITI(1) + PassiveTuning.PreOdor) + [1 PassiveTuning.Odor];
Air2 = Odor(2) + [1  (1+diff(Odor))];
stim_duration = (1+diff(Odor));

% convert Spike rasters to FR
timewindow = 200; % time window over which smoothing occurs - change as needed
t_wid = timewindow;  % width of kernel (in ms)
taxis = -(t_wid*5):(t_wid*5);  % e.g. make a time axis of 1000 ms for a window of 100 ms
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

%% PSTH
PSTH_raw = []; PSTH_zscored = []; PSTH_ratio = []; MeanFRs = [];
% Smooth PSTH - to get a continuous (time-dependent) rate variable
for whichUnit = 1:size(PassiveTuning.RasterOut,1) % each unit
    for whichOdor = 1:size(PassiveTuning.RasterOut,2) % each odor
        for whichLocation = 1:size(PassiveTuning.RasterOut,3) % each location
            for whichRep = 1:size(PassiveTuning.RasterOut,5) % each repeat
                tempPSTH = squeeze(PassiveTuning.RasterOut(whichUnit,whichOdor,whichLocation,:,whichRep))'; % for one cluster, all times for one type
                tempPSTH(:,(Air2(2)+1):end) = [];
                tempPSTH(:,1:(Air1(1)-1)) = [];
                zs = 1000*conv(tempPSTH,gauss_kernel,'same'); %in ms, convert firing rate to Hz (1 ms = 1000 Hz)
                PSTH_raw(whichUnit,:,whichOdor,whichLocation,whichRep) = zs; 
                MeanFRs(whichUnit,1,whichOdor,whichLocation,whichRep) = mean(zs(1,stim_duration+[1:stim_duration])) ; %- ...
                    mean(zs(1,[1:stim_duration]));
                
                % ratio PSTH
                FRo = mean(zs(1,1:stim_duration)); 
                PSTH_ratio(whichUnit,:,whichOdor,whichLocation,whichRep) = (zs - FRo)/FRo;
                
                % z-score
                mu      = mean(zs(1:stim_duration));
                sigma   = std(zs(1:stim_duration));
                PSTH_zscored(whichUnit,:,whichOdor,whichLocation,whichRep) = (zs - mu)/sigma;
            end
        end
    end
end

%% get a sorting order
Unit_attributes = [];
% response strength for center location for each odor
for whichUnit = 1:size(PassiveTuning.RasterOut,1) % each unit
    whichLocation = find(PassiveTuning.Locations==0);
    for m = 1:1:size(PassiveTuning.RasterOut,2) % each odor
        Unit_attributes(whichUnit,m) = mean(MeanFRs(whichUnit,1,m,whichLocation,:));
    end
end

Unit_attributes(:,end+1) = 1:whichUnit;

%% plotting heatmaps
whichOdor = 3; 
Trace2Use = 1; % 1 = RawTraces 2 = RatioTraces 3 = ZscoredTraces
switch Trace2Use
    case 1
        WhichTraces = PSTH_raw;
        range = [0 35];     
    case 2
        WhichTraces = PSTH_ratio;
        range = [-0.1 1];     
    case 3
        WhichTraces = PSTH_zscored;
        range = [-5 70];        
end

% Order by responses of a given Odor, and split left and right bulbs
[~,UnitOrder] = sortrows(Unit_attributes,whichOdor,'ascend'); % order by response strength to center location
                           
% heatmap : all odors, chosen location, average across reps    
whichLocation = find(PassiveTuning.Locations==0);
nStim = numel(PassiveTuning.Odors);
figure('Name',['All Odors, Location = ',num2str(PassiveTuning.Locations(whichLocation))],'NumberTitle','off');
for i = 1:nStim 
    subplot(1,nStim,i); 
    imagesc(mean(WhichTraces(UnitOrder,:,i,whichLocation,:),5), range); 
    set(gca, 'YTick', [], 'XTick', [stim_duration 2*stim_duration], 'XTickLabels', ({}), 'TickDir', 'out');
    colormap(brewermap([],'*RdBu'));
end

% heatmap : given odor, given location, all repeats     
nReps = size(WhichTraces,5);
whichLocation = find(PassiveTuning.Locations==-30);
figure('Name',['All Repeats, Odor = ',num2str(PassiveTuning.Odors(whichOdor)), ', Location = ',num2str(PassiveTuning.Locations(whichLocation))],'NumberTitle','off');
for i = 1:nReps 
    subplot(1,nReps,i); 
    imagesc(WhichTraces(UnitOrder,:,whichOdor,whichLocation,i), range); 
    set(gca, 'YTick', [], 'XTick', [stim_duration 2*stim_duration], 'XTickLabels', ({}), 'TickDir', 'out');
    colormap(brewermap([],'*RdBu'));
end

% heatmap : given odor, all locations, averaged across repeats                 
figure('Name',['All Locations, Odor = ',num2str(PassiveTuning.Odors(whichOdor))],'NumberTitle','off');
nLoc = numel(PassiveTuning.Locations);
for i = 1:nLoc 
    subplot(1,nLoc,i); 
    imagesc(mean(WhichTraces(UnitOrder,:,whichOdor,i,:),5), range); 
    set(gca, 'YTick', [], 'XTick', [stim_duration 2*stim_duration], 'XTickLabels', ({}), 'TickDir', 'out');
    colormap(brewermap([],'*RdBu'));
end
