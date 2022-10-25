%% load GlomSession into work space
%load('/mnt/data/OBImaging/BlackCamK/12-Oct-2022_1/AllGloms.mat');
load('/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/OB imaging/CAMKII_Black/12-Oct-2022_1/AllGloms.mat');

% get stimulus period
preStim = GlomSession.TrialSequence(1,5);
stimulus = (GlomSession.TrialSequence(1,7) - GlomSession.TrialSequence(1,6) + 1);
% endStim = preStim + stimulus;
% postStim = 160;

Frames2Use = [(preStim - stimulus + 1):(preStim + 2*stimulus)];

ZscoredTraces = []; dFTraces = [];
% process all traces - first dF/F then z-score, 
% also reshape - Gloms x time x odors x locations x repeats
for m = 1:length(GlomSession.Odors)
    for n = 1:length(GlomSession.Locations)
        
        whichTrials = find((GlomSession.TrialSequence(:,2) == GlomSession.Odors(m)) & ...
            (GlomSession.TrialSequence(:,1) == GlomSession.Locations(n)));
        
        for trial = 1:numel(whichTrials)
            thisTrial = whichTrials(trial);
            
            for ROI = 1:length(GlomSession.ROI_index)
                
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

whichOdor = 1; 
whichLocation = find(GlomSession.Locations==0);

% Order by responses of a given Odor, and split left and right bulbs
[~,tempOrder] = sortrows(ROI_attributes,whichOdor+1,'descend'); % order by response strength
GlomsOrdered = ROI_attributes(tempOrder,:);
% split left and right into chunks, with reversed orders
GlomOrder = tempOrder(vertcat( find(GlomsOrdered(:,1)<MidLine),...
                               flipud(find(GlomsOrdered(:,1)>MidLine)) ));

%% heatmap: all odors, center location, averaged across reps
range = [-5 80];
nStim = numel(GlomSession.Odors);
for i = 1:nStim
    subplot(1,nStim,i); 
    imagesc(mean(ZscoredTraces(GlomOrder,:,i,whichLocation,:),5), range); 
    set(gca, 'YTick', [], 'XTick', []);
    colormap(brewermap([],'*RdBu'));
end
                           
%% heatmap : given odor, given location, all repeats

for whichOdor = 2%1:nStim
    
    % get sorting order
    
    % Order by responses of a given Odor, and split left and right bulbs
    [~,tempOrder] = sortrows(ROI_attributes,whichOdor+1,'descend'); % order by response strength
    GlomsOrdered = ROI_attributes(tempOrder,:);
    % split left and right into chunks, with reversed orders
    GlomOrder = tempOrder(vertcat( find(GlomsOrdered(:,1)<MidLine),...
        flipud(find(GlomsOrdered(:,1)>MidLine)) ));
    
    range = [-5 80];
    nReps = size(ZscoredTraces,5);
    figure;
    
    whichLocation = find(GlomSession.Locations==-10);
    
    % plot the average
    subplot(1,nReps+1,1);
    imagesc(mean(ZscoredTraces(GlomOrder,:,whichOdor,whichLocation,:),5), range);
    set(gca, 'YTick', [], 'XTick', []);
    colormap(brewermap([],'*RdBu'));
    
    for i = 1:nReps
        subplot(1,nReps+1,i+1);
        imagesc(ZscoredTraces(GlomOrder,:,whichOdor,whichLocation,i), range);
        set(gca, 'YTick', [], 'XTick', []);
        colormap(brewermap([],'*RdBu'));
    end

end

%% heatmap : given odor, all locations, averaged across repeats        
for whichOdor = 1:nStim
    
     % get sorting order
    whichLocation = find(GlomSession.Locations==0);
    
    % Order by responses of a given Odor, and split left and right bulbs
    [~,tempOrder] = sortrows(ROI_attributes,whichOdor+1,'descend'); % order by response strength
    GlomsOrdered = ROI_attributes(tempOrder,:);
    % split left and right into chunks, with reversed orders
    GlomOrder = tempOrder(vertcat( find(GlomsOrdered(:,1)<MidLine),...
        flipud(find(GlomsOrdered(:,1)>MidLine)) ));
    
    figure;
    nLoc = numel(GlomSession.Locations);
    for i = 1:nLoc
        subplot(1,nLoc,i);
        imagesc(mean(ZscoredTraces(GlomOrder,:,whichOdor,i,:),5), range);
        set(gca, 'YTick', [], 'XTick', []);
        colormap(brewermap([],'*RdBu'));
    end
end