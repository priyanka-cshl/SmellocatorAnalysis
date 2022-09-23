function [Ratio, validTrials] = GetRatioImages(FolderPath)

if ~exist('FolderPath','var')
    %FolderPath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/OB imaging/26-Jul-2022_4';
    FolderPath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/OB imaging/04-Aug-2022';
    
end

% Load all trials and construct Average Ratio Images for each
nTrials = numel(dir([FolderPath,filesep,'frameTimes*.mat']));
threshold = 1500; % for digitizing triggers from the analog file
validTrials = zeros(nTrials,2);
Ratio = zeros(540,640,nTrials);

for i = 1:nTrials
    disp(i);
    % get info about this trial
    % get frame times as saved by Widefield Imager
    load([FolderPath,filesep,'frameTimes_',num2str(i, '%04i'),'.mat'],'preStim','imgSize'); %, 'imgSize','frameTimes');
    
    % read frame/stimulus triggers from analog file to get StimON period
    [~,Analog] = loadRawData(fullfile(FolderPath,['Analog_',num2str(i),'.dat']),'Analog');
    Analog = double(Analog);
    FrameStarts = sort([find(diff(Analog(2,:))>threshold) find(diff(Analog(3,:))>threshold)]); % ignore blue vs violet channels
    StimStart = find(diff(Analog(4,:))>threshold);
    StimStop = find(diff(Analog(4,:))<-threshold);
    validTrials(i,2) = mode(Analog(5,StimStart:StimStop));
    
    if ~isempty(StimStart)
        StimDuration = (find(FrameStarts>StimStop,1,'first') - find(FrameStarts>StimStart,1,'first'));
    end
    
    if preStim > 10 % valid trial
        validTrials(i,1) = 1;
        
        % load Imaging Data
        tag = num2str(i);
        appender = [];
        for x = 1:(4-length(tag))
            appender = [appender,'0'];
        end
        tag = ['Frames*',appender,tag,'.*'];
        MyFile = dir([FolderPath,filesep,tag]);
        [~, MyStack] = loadRawData(fullfile(FolderPath,MyFile.name),'frames',[],imgSize);
        
        Baseline = mean(MyStack(:,:,:,(preStim-10):preStim),4);
        Stimulus = mean(MyStack(:,:,:,preStim+(1:StimDuration)),4);
        %Stimulus = mean(MyStack(:,:,:,preStim+(11:StimDuration)),4);
        Ratio(:,:,i) = (Stimulus - Baseline)./Baseline;
    end
    
end

end