function [GlomTraces, ROIList] = GetGlomTraces(ROIImage, ROIList, FolderPath, imgSize, TrialSettings, FrameWindow)

if nargin<6
    FrameWindow = [10 10];
end

if isempty(ROIList)
    ROIList = unique(ROIImage);
    ROIList(ROIList==0,:) = [];
end

GlomTraces = [];
for thisTrial = 1:size(TrialSettings,1)
    
    TrialNum        = TrialSettings(thisTrial,1);
    PreStim         = TrialSettings(thisTrial,2);
        
    if FrameWindow(1)>0
        BaselineFrames = (PreStim - FrameWindow(1) + 1):PreStim;
    else
        BaselineFrames = 1:PreStim;
    end
    
    % load Imaging Data
    tag = num2str(TrialNum);
    appender = [];
    for x = 1:(4-length(tag))
        appender = [appender,'0'];
    end
    tag = ['Frames*',appender,tag,'.*'];
    
    MyFile = dir([FolderPath,filesep,tag]);
    [~, MyStack] = loadRawData(fullfile(FolderPath,MyFile.name),'frames',[],imgSize);
    nPix = numel(MyStack(:,:,:,1));
    MyStack = reshape(MyStack,nPix,[]); % linearized stack
    for j = 1:numel(ROIList)
        ROImask = find(ROIImage==ROIList(j));
        thisGlomTrace = mean(MyStack(ROImask,:));
        GlomTraces(j,:,1) = thisGlomTrace;
        Fo = mean(thisGlomTrace(BaselineFrames));
        GlomTraces(j,:,2) = (thisGlomTrace - Fo)/Fo;
    end
    
end

