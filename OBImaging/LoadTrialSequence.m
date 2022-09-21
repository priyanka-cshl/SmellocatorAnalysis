function [MergedTrialList, Locations, Odors, Reps, ImageSize] = LoadTrialSequence(BehaviorPath,ImagingPath)

% trial list from the behavior GUI
load(BehaviorPath);
TrialList_Behavior = session_data.TrialSequence;

% trial list from the imaging side
nTrials = numel(dir([ImagingPath,filesep,'frametimes*.mat']));
threshold = 1500; % for digitizing triggers from the analog file
TrialList_Imaging = zeros(nTrials,6);
ImageSize = 0;
for i = 1:nTrials
    % get info about this trial
    TrialList_Imaging(i,1) = i;
    
    % get frame times as saved by Widefield Imager
    load([ImagingPath,filesep,'frameTimes_',num2str(i, '%04i'),'.mat'],'preStim','postStim','imgSize'); %, 'imgSize','frameTimes');
    
    % read frame/stimulus triggers from analog file to get StimON period
    [~,Analog] = loadRawData(fullfile(ImagingPath,['Analog_',num2str(i),'.dat']),'Analog');
    Analog = double(Analog);
    FrameStarts = sort([find(diff(Analog(2,:))>threshold) find(diff(Analog(3,:))>threshold)]); % ignore blue vs violet channels
    StimStart = find(diff(Analog(4,:))>threshold);
    StimStop = find(diff(Analog(4,:))<-threshold);
    TrialList_Imaging(i,2) = mode(Analog(5,StimStart:StimStop));
    
    if ~isempty(StimStart)
        StimDuration = (find(FrameStarts>StimStop,1,'first') - find(FrameStarts>StimStart,1,'first'));
    end
    
    TrialList_Imaging(i,3:6) = [preStim find(FrameStarts>StimStart,1,'first') (find(FrameStarts>StimStop,1,'first')-1) postStim];
    
    if i == 3
        ImageSize = imgSize;
    end
end

% Match the two trial lists
% Assume first two trials in the Imaging side are weird
% check by fitting a line to motor locations from behavior side and imaging side (analog data)
ExtraTrials = size(TrialList_Imaging,1) - size(TrialList_Behavior,1);

if ExtraTrials
    [~,gof] = fit(TrialList_Imaging((ExtraTrials+1):end,2),TrialList_Behavior(:,1),'poly1');
    if gof.rsquare<0.9
        error('Unable to match Imaging and Behavior Trial Sequence');
    else
        MergedTrialList = [TrialList_Behavior TrialList_Imaging(ExtraTrials+1:end,:)];
    end
end

Locations = unique(TrialList_Behavior(:,1));
Odors     = unique(TrialList_Behavior(:,2));
Reps      = size(TrialList_Behavior,1)/numel(Odors)/numel(Locations);

end




    

