function [GlomSession] = GetAllGlomTraces(FolderPath)

clear GlomSession

% Get Trial Sequence
BehaviorFile = dir([FolderPath,'*_o*.mat']);
[TrialSequence, GlomSession.Locations, GlomSession.Odors, GlomSession.Reps, ImageSize] = ...
    LoadTrialSequence(fullfile(FolderPath,BehaviorFile.name), FolderPath);

% Get all ROIs
clear ROIs
if exist(fullfile(FolderPath,'GlomerularMasks.mat'))
    load(fullfile(FolderPath,'GlomerularMasks.mat'), 'ROIs');
end
ROIImage = ROIs;
clear ROIs
ROIList = unique(ROIImage);
ROIList(ROIList==0,:) = [];

% get traces for each glomerulus, each trial
GlomTraces = [];
for thisTrial = 1:size(TrialSequence,1) % every trial
        
    TrialNum        = TrialSequence(thisTrial,3);
    PreStim         = TrialSequence(thisTrial,5);
        
    % load Imaging Data
    tag = num2str(TrialNum);
    appender = [];
    for x = 1:(4-length(tag))
        appender = [appender,'0'];
    end
    tag = ['Frames*',appender,tag,'.*'];
    
    MyFile = dir([FolderPath,filesep,tag]);
    [~, MyStack] = loadRawData(fullfile(FolderPath,MyFile.name),'frames',[],ImageSize);
    nPix = numel(MyStack(:,:,:,1));
    MyStack = reshape(MyStack,nPix,[]); % linearized stack
    
    for j = 1:numel(ROIList) % every ROI
        ROImask = find(ROIImage==ROIList(j));
        thisGlomTrace = mean(MyStack(ROImask,:));
        GlomTraces(j,1:length(thisGlomTrace),thisTrial) = thisGlomTrace;
    end
end

% All outputs
GlomSession.Traces = GlomTraces;
GlomSession.ROI_index = ROIList;
GlomSession.TrialSequence = TrialSequence;
    

