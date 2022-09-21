function [Ratio] = GetRatioImage(FolderPath, imgSize, TrialSettings, FrameWindow)

if nargin<3
    FrameWindow = [10 10];
end

for thisTrial = 1:size(TrialSettings,1)
    
    TrialNum        = TrialSettings(thisTrial,1);
    PreStim         = TrialSettings(thisTrial,2);
    StimStop        = PreStim + (TrialSettings(thisTrial,4) - TrialSettings(thisTrial,3) + 1);
    
    if FrameWindow(1)>0
        BaselineFrames = (PreStim - FrameWindow(1) + 1):PreStim;
    else
        BaselineFrames = 1:PreStim;
    end
    
    if FrameWindow(2)>0
        StimulusFrames = (StimStop - FrameWindow(2) + 1):StimStop;
    else
        StimulusFrames = (PreStim + 1):StimStop;
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
    
    Baseline = mean(MyStack(:,:,:,BaselineFrames),4);
    Stimulus = mean(MyStack(:,:,:,StimulusFrames),4);
    Ratio(:,:,thisTrial) = (Stimulus - Baseline)./Baseline;
    
end

