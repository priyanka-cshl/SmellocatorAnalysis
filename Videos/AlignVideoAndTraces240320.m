function [handles] = AlignVideoAndTraces240319(slowby, write_video)

if nargin<1
    slowby = 1;
    write_video = 0;
end
if nargin<2
    write_video = 0;
end

passive_replay = 1;

%% FilePaths
RootVideo = '/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/Movies/Raw';
VideoFolder{1} = fullfile(RootVideo,'20181025_cam1_N8');
VideoFolder{2} = fullfile(RootVideo,'20181025_cam2_N8');
filetag{1} = 'fc3_save_2018-10-24-175140-'; %0000';
filetag{2} = 'fc3_save_2018-10-24-175139-'; %0000';
DataFile = fullfile('/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/N8/N8_20181024_r1.mat');
%slowby = 1;

%% Read the Data File
SetupSmellocatorGlobals;
[MyData, MySettings, DataTags] = ReadSessionData(DataFile);
[Trials,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);
[AllTFs] = SmellocatorTransferFunctions;
AllTFs = -AllTFs;
MaxMotorLocation = 100;

%% find the time periods when camera was acquiring in sync
camera_start = MySettings(find(diff(MySettings(1:end,3))==1)+1,1); % +1 because of diff
camera_stop = MySettings(find(diff(MySettings(1:end,3))==-1)+1,1); % +1 because of diff
% pick period of interest, in case there are multiple video runs
camera_start(1,:) = []; camera_stop(1,:) = [];
% ignore camera triggers outside the period of interest
CameraCols = find(cellfun(@isempty,regexp(DataTags,'Pgrey*'))==0)';
MyData(find(MyData(:,1)<camera_start),CameraCols) = 0*MyData(find(MyData(:,1)<camera_start),CameraCols);
MyData(find(MyData(:,1)>camera_stop),CameraCols) = 0*MyData(find(MyData(:,1)>camera_stop),CameraCols);
% get frame time stamps
FrameTimeStamps_Cam1 = MyData(find(diff(MyData(:,CameraCols(1)))==1),1);
FrameTimeStamps_Cam2 = MyData(find(diff(MyData(:,CameraCols(2)))==1),1);

% define stretch of time you want to use for the video (in seconds)
tstart = 297; tstop = 311; % tstart = 291; tstop = 319;
%tstart = 304; tstop = 311; % tstart = 291; tstop = 319;
%tstart = 312; tstop = 320; % tstart = 291; tstop = 319;
% tstart = 274; tstop = 290; % tstart = 291; tstop = 319;
%tstart = 290; tstop = 305; % tstart = 291; tstop = 319;
FrameStartIdx_Cam1 = find(FrameTimeStamps_Cam1>=tstart,1);
FrameStartIdx_Cam2 = find(FrameTimeStamps_Cam2>=tstart,1);
FrameStopIdx_Cam1 = find(FrameTimeStamps_Cam1<=tstop,1,'last');
FrameStopIdx_Cam2 = find(FrameTimeStamps_Cam2<=tstop,1,'last');
frames_to_show = FrameStopIdx_Cam1 - FrameStartIdx_Cam1 + 1;
framerate = round(frames_to_show/(tstop-tstart));
timewindow = [-0.25 (tstop-tstart+1)]; %[-0.25 2]; % seconds
%frames_to_show = 100;

%% trialstarts and ends
trialstarts = find(diff(MyData(:,5))>0);
trialends = find(diff(MyData(:,5))<0);
trialstarts(MyData(trialstarts,1)<FrameTimeStamps_Cam1(1),:) = [];
trialstarts(MyData(trialstarts,1)>FrameTimeStamps_Cam1(end),:) = [];
trialends(trialends<trialstarts(1),:) = [];
FrameList = [];
for trials = 1:numel(trialstarts)
    whichframe = find(FrameTimeStamps_Cam1<=MyData(trialstarts(trials),1),1,'last'); % last frame before trial start
    FrameList(trials,1) = whichframe;
    whichframe = find(FrameTimeStamps_Cam1<=MyData(trialends(trials),1),1,'last'); % last frame before trial end
    FrameList(trials,2) = whichframe;
    FrameList(trials,3) = MyData(trialends(trials),4); % Lever Position
end

%% make a list of frames from camera 1 where lever is in trial start state
trialstartframes = [];
for n = 1:numel(FrameTimeStamps_Cam1)
    idx = find(MyData(:,1)==FrameTimeStamps_Cam1(n),1);
    if MyData(idx,4) > 4.9 && MyData(idx,4) < 5
        trialstartframes = [trialstartframes; n];
    end
end
trialstartframes = unique(trialstartframes);
if numel(trialstartframes)<frames_to_show
    trialstartframes = repmat(trialstartframes,ceil(frames_to_show/numel(trialstartframes)),1);
end

%% Video writing related initializations
if write_video
    %writerObj = VideoWriter('test_video','MPEG-4'); % VideoWriter('test_video','avi'); %
    %writerObj = VideoWriter('test_video_slow10.avi');
    %writerObj = VideoWriter(['N8_withPassiveReplay_slow',num2str(slowby),'.avi']);
    writerObj = VideoWriter(['N8_withPassiveReplay_slow',num2str(slowby)],'MPEG-4');
    writerObj.FrameRate = framerate/slowby;
    open(writerObj);
end

%% Video display related initializations
% display parameters
a = 0.8; %0.5; % proportion occupied by the video image - choose between 0 and 1
margin = 0.5;
edge = 0.05;

% load and display one frame
frame = imread(fullfile(VideoFolder{1},[filetag{1},'0000.png']));
frame = imresize(frame,0.5); % half the image size
ratio = size(frame,1)/size(frame,2);

% to get the resolution of the monitor minus the menu elements
video_figure = figure('name','video');
%set(video_figure, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 1]);
set(video_figure, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
myframe = getframe(video_figure);
framedimensions = size(myframe.cdata);
screenratio = framedimensions(2)/framedimensions(1);
close (gcf);
pause(0.5);

%% Actual figure combining video and traces
h = figure('position', [1 1 framedimensions(2) framedimensions(1)], 'color', 'w');
pause(0.5);

%% Behavior data
% initialize the data axes and plot the session data
handles.BehaviorPlot = subplot(2,2,[1 2]);
[timestamps,Lever,Sniffs,Licks,Rewards,Trial,TZ] = PrepBehaviorForPlotting(MyData, DataTags);
PlotBehavior(timestamps,Lever,Sniffs,Licks,Rewards,Trial,TZ,5,'plotTF',0);
set(handles.BehaviorPlot, 'position', [-0.2+margin/2, 0.7, 0.85, 0.2]);
set(handles.BehaviorPlot,'XTick',handles.BehaviorPlot.XLim(1):2:handles.BehaviorPlot.XLim(2),'XTickLabel',num2str([(handles.BehaviorPlot.XLim(1):2:handles.BehaviorPlot.XLim(2))-handles.BehaviorPlot.XLim(1)]'))
set(handles.BehaviorPlot,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
handles.BehaviorPlot.XLim = tstart + timewindow;
handles.BehaviorPlot.FontSize = 12;
handles.BehaviorPlot.FontWeight = 'bold';
set(handles.BehaviorPlot,'YTick',[0 5],'YTickLabel',{'min' 'max'},'TickDir','out','Box','off');
handles.BehaviorPlot.TickLength(1) = 0.005;

%% initialize another axes object to plot the timestamp bar
handles.TimeBar = axes;
handles.TimeBar.Position = handles.BehaviorPlot.Position;
handles.TimeBar.Position(2) = handles.TimeBar.Position(2) - 0.005;
handles.TimeBar.Position(4) = handles.TimeBar.Position(4) + 0.01;
axes(handles.TimeBar);
set(gca, 'Color', 'none');
handles.TimeBar.YLim = handles.BehaviorPlot.YLim;
handles.TimeBar.XLim = tstart + timewindow;
set(handles.TimeBar,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.TimeBar,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual','TickDir','out','Box','off');
handles.bar = line([tstart tstart],handles.TimeBar.YLim,'color',[239 41 41]./256,'LineWidth',2);
handles.TimeBar.FontSize = 14;
handles.TimeBar.TickLength(1) = 0.005;
handles.TimeBar.YColor = 'none';
handles.TimeBar.XColor = 'none';

%% initialize Camera1 axes with one frame
handles.Camera1 = subplot(2,2,3);
axes(handles.Camera1);
handles.Camera1.LineWidth = 2;
handles.frame1 = imagesc(frame,'parent',handles.Camera1);
colormap('gray');
set(handles.Camera1,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.Camera1,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
set(handles.Camera1, 'position', [margin/10, edge , a/2, screenratio*ratio*a/2]);
% add a trial state label
% textposition = handles.Camera1.Position;
% handles.trialstate = text(625,125,'ITI','color','w','Fontsize', 28, 'Fontweight','bold', 'HorizontalAlignment', 'right');

%% initialize Camera2 axes with one frame
handles.Camera2 = subplot(2,2,4);
axes(handles.Camera2);
handles.Camera2.LineWidth = 2;
handles.frame2 = imagesc(frame,'parent',handles.Camera2);
colormap('gray');
set(handles.Camera2,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.Camera2,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
set(handles.Camera2, 'position', [(screenratio*ratio*a/2)-margin/7, edge , a/2, screenratio*ratio*a/2]);

%% add another axes to mark the targetzone boundaries
zonehalfwidth = 0.0625*handles.Camera1.Position(3);
zoneheight = 0.06;
handles.TZLims = axes;
handles.TZLims.Position(1) = handles.Camera1.Position(1) + (handles.Camera1.Position(3))/2 - zonehalfwidth;
handles.TZLims.Position(2) = handles.Camera1.Position(2) + (handles.Camera1.Position(4));
handles.TZLims.Position(3) = 2*zonehalfwidth;
handles.TZLims.Position(4) = zoneheight/2;
set(handles.TZLims,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.TZLims,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual','Box','on');
handles.TZLims.Color = [1 1 1];

if passive_replay
%% add another axes to mark the targetzone boundaries on the passive replay movie
handles.FakeTZLims = axes;
handles.FakeTZLims.Position(1) = handles.Camera2.Position(1) + (handles.Camera2.Position(3))/2 - zonehalfwidth;
handles.FakeTZLims.Position(2) = handles.Camera2.Position(2) + (handles.Camera2.Position(4));
handles.FakeTZLims.Position(3) = 2*zonehalfwidth;
handles.FakeTZLims.Position(4) = zoneheight/2;
set(handles.FakeTZLims,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.FakeTZLims,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual','Box','on');
handles.FakeTZLims.Color = [1 1 1];
end

%% some housekeeping
frame_offset = 0;
frameAID = FrameStartIdx_Cam1 - 1;
frameBID = 1;

% rectangle params to capture just the 2 cameras
buffermargin = 0.01;
leftlim     = handles.Camera1.Position(1) - buffermargin;
rightlim    = handles.Camera2.Position(1) + handles.Camera2.Position(3) + buffermargin;
bottomlim   = handles.Camera1.Position(2) - buffermargin;
toplim      = handles.TZLims.Position(2) + handles.TZLims.Position(4) + buffermargin;
R           = [leftlim bottomlim (rightlim - leftlim) (toplim - bottomlim)];
figsize = get(gcf,'Position');
Rpix        = round(R.*([figsize(3) figsize(4) figsize(3) figsize(4)]));

trialcount = 0;
trialphase = -1; % ITI

ManifoldColors(1,:) = [0.8 0.8 0.8 0.6]; % ITI
ManifoldColors(2,:) = [1 1 1 0.6]; %Trial
ManifoldColors(3,:) = [0.9490 0.9176 0.6980 0.6]; % TZ
ManifoldColors(2,:) = [1 1 1 0.0]; %Trial
ManifoldColors(3,:) = [1 1 1 0.0]; %Trial
%handles.OdorBox.FaceColor = ManifoldColors(1,:);

%% actual display of frames
for i = 1:frames_to_show
    frameAID = frameAID + 1;
    %frameBID = frameBID + 1;
    
    %% update frames
    %axes(handles.Camera1);
    frameA = imread(fullfile(VideoFolder{1},[filetag{1},sprintf('%04d',frameAID),'.png']));
    handles.frame1.CData = imresize(frameA,0.5);
    
    %axes(handles.Camera2);
    if ~passive_replay
        frameB = imread(fullfile(VideoFolder{2},[filetag{2},num2str(frameBID),'.png']));
    else
        %frameB = imread(fullfile(VideoFolder{1},[filetag{1},sprintf('%04d',trialstartframes(frameBID)),'.png']));
        frameB = imread(fullfile(VideoFolder{1},[filetag{1},sprintf('%04d',4646),'.png']));
        % replace top parts with frame A
        frameB(1:200,:) = frameA(1:200,:);
        handles.frame2.CData = imresize(frameB,0.5);
    end
    
    %% update time bar position
    tstart = FrameTimeStamps_Cam1(frameAID);
    handles.bar.XData = [tstart tstart];
    
    %% update trial state string and TZLims color
    % for updating trial states etc
    % find the behavioral sample index given current timestamp
    idx = find(abs(MyData(:,1)-tstart)==min(abs(MyData(:,1)-tstart)),1);
    
    % detect trial start (trial col = 5)
    if trialphase < 1 && MyData(idx,5) > 0
        trialcount = trialcount + 1;
        handles.trialstate.String = ['Trial ', num2str(trialcount)];
        trialphase = 1;
        handles.TZLims.Color = [0.8 0.8 0.8]; % grey
        if passive_replay
            handles.FakeTZLims.Color = [0.8 0.8 0.8];
        end
        %handles.trialstate.Color = [1 1 1];
        %handles.OdorBox.FaceColor = ManifoldColors(2,:);
        %handles.OdorBox.FaceColor = [0.9490 0.9176 0.6980 0.8];
        %FrameList(trialcount,1) = frameAID;
    end
        
    % detect trial end (trial col = 5)
    if trialphase > 0 && MyData(idx,5) == 0
        handles.trialstate.String = 'ITI';
        trialphase = -1;
        handles.TZLims.Color = [1 1 1];
        if passive_replay
            handles.FakeTZLims.Color = [1 1 1];
        end
        %handles.trialstate.Color = [1 1 1];
        %handles.OdorBox.FaceColor(4) = 0.5;
        %handles.OdorBox.FaceColor = ManifoldColors(1,:);
        %FrameList(trialcount,2) = frameAID;
    end
    
    % detect trial triggers
    if trialphase == -1 % in ITI
        if MyData(idx,4) >= 4.8
            %handles.trialstate.String = 'Trigger';
            trialphase = 0;
        end
    end
    
    %% update TZlims box based on whether the lever is in target or not
    if  trialphase>0 % trial is ON
        if MyData(idx,6)>0 && mode(MyData(idx-19:idx,7))>0 % in trial, in TargetZone
%             handles.TZLims.Color = [1 1 0];
%             handles.trialstate.Color = [1 1 0];
            handles.TZLims.Color = [242 234 178]./255;
            if passive_replay
                handles.FakeTZLims.Color = [242 234 178]./255;
            end
            %handles.trialstate.Color = [242 234 178]./255;
            %handles.OdorBox.FaceColor = ManifoldColors(3,:);
        else % in trial, not in targetzone
            handles.TZLims.Color = [0.8 0.8 0.8];
            if passive_replay
                handles.FakeTZLims.Color = [0.8 0.8 0.8];
            end
            %handles.trialstate.Color = [1 1 1];
            %handles.OdorBox.FaceColor = ManifoldColors(2,:);
        end
    end
        
    if write_video
        myframe = getframe(h,Rpix);
        %for foo = 1:slowby
            writeVideo(writerObj,myframe);
        %end
    else
        pause(1/(framerate/slowby));
    end
    
end

if write_video
    close(writerObj);
    delete(writerObj);
end

close(h);
%%
% figure('position', [1 1 framedimensions(2) framedimensions(1)], 'color', 'w');
% colormap('gray');
% for f = 1:size(FrameList,1)-1
%     frameID = FrameList(f,1);
%     frame = imread(fullfile(VideoFolder{1},[filetag{1},sprintf('%04d',frameID),'.png']));
%     frame = imresize(frame,0.5);
%     subplot(size(FrameList,1)-1,2,(f*2)-1)
%     imagesc(frame);
%     set(gca,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
%     set(gca,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
% 
%     frameID = FrameList(f,2);
%     frame = imread(fullfile(VideoFolder{1},[filetag{1},sprintf('%04d',frameID),'.png']));
%     frame = imresize(frame,0.5);
%     subplot(size(FrameList,1)-1,2,(f*2))
%     imagesc(frame);
%     set(gca,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
%     set(gca,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
% 
% end

end



