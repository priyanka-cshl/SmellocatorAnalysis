%% load GlomSession into work space
load('/mnt/data/OBImaging/BlackCamK/14-Oct-2022_1/AllGloms.mat');
%load('/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/OB imaging/CAMKII_Black/14-Oct-2022_1/AllGloms.mat');

% get stimulus period
preStim = GlomSession.TrialSequence(1,5);
stimulus = (GlomSession.TrialSequence(1,7) - GlomSession.TrialSequence(1,6) + 1);
% endStim = preStim + stimulus;
% postStim = 160;

% which frames to keep
Frames2Use = [(preStim - stimulus + 1):(preStim + 2*stimulus)];

RawTraces = []; ZscoredTraces = []; dFTraces = [];
% process all traces - first dF/F then z-score, 
% also reshape - Gloms x time x odors x locations x repeats
for m = 1:length(GlomSession.Odors)
    for n = 1:length(GlomSession.Locations)
        
        whichTrials = find((GlomSession.TrialSequence(:,2) == GlomSession.Odors(m)) & ...
            (GlomSession.TrialSequence(:,1) == GlomSession.Locations(n)));
        
        for trial = 1:numel(whichTrials)
            thisTrial = whichTrials(trial);
            
            for ROI = 1:length(GlomSession.ROI_index)
                
                % Raw
                RawTraces(ROI,:,m,n,trial) = GlomSession.Traces(ROI, Frames2Use, thisTrial);
                
                % dF/F
                Fo = mean(GlomSession.Traces(ROI,1:(GlomSession.TrialSequence(thisTrial,5)), thisTrial));
                dFTraces(ROI,:,m,n,trial) = (GlomSession.Traces(ROI, Frames2Use, thisTrial) - Fo)/Fo;
                
                % z-score
                mu = mean(GlomSession.Traces(ROI,1:(GlomSession.TrialSequence(thisTrial,5)), thisTrial));
                sigma = std(GlomSession.Traces(ROI,1:(GlomSession.TrialSequence(thisTrial,5)), thisTrial));
                
                ZscoredTraces(ROI,:,m,n,trial) = (GlomSession.Traces(ROI, Frames2Use, thisTrial) - mu)/sigma;
                
            end
        end
    end
end

%% get a glomerular order
% by location left vs. right bulb
% and response strength for center location for each odor
ROI_attributes = [];
MidLine = 308;
for ROI = 1:length(GlomSession.ROI_index)
    % find approx centroid on the x-axis
    ROI_attributes(ROI,1) = mean(ceil(find(GlomSession.ROImasks==GlomSession.ROI_index(ROI)))/size(GlomSession.ROImasks,1));
    % response strengths for each odor 
    Frames2Use = stimulus+(1:stimulus); 
    whichLocation = find(GlomSession.Locations==0);
    for m = 1:length(GlomSession.Odors)
        ROI_attributes(ROI,1+m) = mean(mean(ZscoredTraces(ROI,:,m,whichLocation,:)));
    end
end

ROI_attributes(:,end+1) = 1:ROI;

        
%% plotting heatmaps
whichOdor = 3; 
Trace2Use = 1; % 1 = RawTraces 2 = dFTraces 3 = ZscoredTraces
switch Trace2Use
    case 1
        WhichTraces = RawTraces;
        range = [0 15000];     
    case 2
        WhichTraces = dFTraces;
        range = [-0.1 1];     
    case 3
        WhichTraces = ZscoredTraces;
        range = [-5 70];        
end

% Order by responses of a given Odor, and split left and right bulbs
[~,tempOrder] = sortrows(ROI_attributes,whichOdor+1,'descend'); % order by response strength to center location
GlomsOrdered = ROI_attributes(tempOrder,:);
% split left and right into chunks, with reversed orders
GlomOrder = tempOrder(vertcat( find(GlomsOrdered(:,1)<MidLine),...
                               flipud(find(GlomsOrdered(:,1)>MidLine)) ));

% heatmap : all odors, chosen location, average across reps    
whichLocation = find(GlomSession.Locations==0);
nStim = numel(GlomSession.Odors);
figure;
for i = 1:nStim 
    subplot(1,nStim,i); 
    imagesc(mean(WhichTraces(GlomOrder,:,i,whichLocation,:),5), range); 
    set(gca, 'YTick', [], 'XTick', []);
    colormap(brewermap([],'*RdBu'));
end

% heatmap : given odor, given location, all repeats     
nReps = size(WhichTraces,5);
whichLocation = find(GlomSession.Locations==0);
figure;
for i = 1:nReps 
    subplot(1,nReps,i); 
    imagesc(WhichTraces(GlomOrder,:,whichOdor,whichLocation,i), range); 
    set(gca, 'YTick', [], 'XTick', []);
    colormap(brewermap([],'*RdBu'));
end

% heatmap : given odor, all locations, averaged across repeats                 
figure;
nLoc = numel(GlomSession.Locations);
for i = 1:nLoc 
    subplot(1,nLoc,i); 
    imagesc(mean(WhichTraces(GlomOrder,:,whichOdor,i,:),5), range); 
    set(gca, 'YTick', [], 'XTick', []);
    colormap(brewermap([],'*RdBu'));
end

