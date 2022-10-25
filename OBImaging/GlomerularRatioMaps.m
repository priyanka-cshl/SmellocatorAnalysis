%% load GlomSession into work space
%load('/mnt/data/OBImaging/BlackCamK/12-Oct-2022_1/AllGloms.mat');
load('/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/OB imaging/CAMKII_Black/14-Oct-2022_1/AllGloms.mat');

% get stimulus period
preStim = GlomSession.TrialSequence(1,5);
stimulus = (GlomSession.TrialSequence(1,7) - GlomSession.TrialSequence(1,6) + 1);

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
        ROI_attributes(ROI,1+m) = max(mean(ZscoredTraces(ROI,:,m,whichLocation,:),5));
    end
end

ROI_attributes(:,end+1) = 1:ROI;

%% Get a Ratio Image
RatioImage = [];
for odor = 1:3
    allROIs = GlomSession.ROImasks;
    for roi = 1:numel(GlomSession.ROI_index)
        whichPixels = find(allROIs == roi);
        
        allROIs(whichPixels) = ROI_attributes(roi,2+odor);
    end
    RatioImage(:,:,odor) = allROIs;
end

%%
for odor = 1:3
    figure;
    MyImage = RatioImage(:,:,odor);
    MyImage(MyImage==0) = - 100;
    imagesc(MyImage,[-100 150]); 
    switch odor
        case 1
            colormap(brewermap([],'Reds'));
        case 2
            colormap(brewermap([],'Blues'));
        case 3
            colormap(brewermap([],'Greens'));
    end
end
