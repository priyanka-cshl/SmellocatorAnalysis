function [handles] = AlignVideoAndTraces240318(write_video)

if nargin<1
    write_video = 0;
end

%% FilePaths
RootVideo = '/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/Movies/Raw';
VideoFolder{1} = fullfile(RootVideo,'20181025_cam1_N8');
VideoFolder{2} = fullfile(RootVideo,'20181025_cam2_N8');
filetag{1} = 'fc3_save_2018-10-24-175140-'; %0000';
filetag{2} = 'fc3_save_2018-10-24-175139-'; %0000';
DataFile = fullfile('/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/N8/N8_20181024_r1.mat');
slowby = 2;

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
FrameStartIdx_Cam1 = find(FrameTimeStamps_Cam1>=tstart,1);
FrameStartIdx_Cam2 = find(FrameTimeStamps_Cam2>=tstart,1);
FrameStopIdx_Cam1 = find(FrameTimeStamps_Cam1<=tstop,1,'last');
FrameStopIdx_Cam2 = find(FrameTimeStamps_Cam2<=tstop,1,'last');
frames_to_show = FrameStopIdx_Cam1 - FrameStartIdx_Cam1 + 1;
framerate = frames_to_show/(tstop-tstart+1);
timewindow = [-0.25 (tstop-tstart+1)]; %[-0.25 2]; % seconds

%% Video writing related initializations
if write_video
    writerObj = VideoWriter('test_video','MPEG-4');
    writerObj.FrameRate = framerate;
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

% initialize the data axes and plot the session data
handles.H2 = subplot(2,2,[1 2]);
[timestamps,Lever,Sniffs,Licks,Rewards,Trial,TZ] = PrepBehaviorForPlotting(MyData, DataTags);
PlotBehavior(timestamps,Lever,Sniffs,Licks,Rewards,Trial,TZ,5,'plotTF',0);

% initialize video axes with one frame
handles.H1 = subplot(2,2,3);
axes(handles.H1);
handles.H1.LineWidth = 2;
handles.frame = imagesc(frame,'parent',handles.H1);
colormap('gray');
set(handles.H1,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.H1,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
%set(handles.H1, 'position', [margin/2, edge*2 , a, a*ratio*framedimensions(2)/framedimensions(1)]); % left, bottom, width, height relative to the bottom left corner
set(handles.H1, 'position', [margin/10, edge , a/2, screenratio*ratio*a/2]);

% initialize video axes with one frame
handles.H1B = subplot(2,2,4);
axes(handles.H1B);
handles.H1B.LineWidth = 2;
handles.frameB = imagesc(frame,'parent',handles.H1B);
colormap('gray');
set(handles.H1B,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.H1B,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
set(handles.H1B, 'position', [(screenratio*ratio*a/2)-margin/10, edge , a/2, screenratio*ratio*a/2]);

set(handles.H2, 'position', [-0.2+margin/2, 0.7, 0.85, 0.2]);
set(handles.H2,'XTick',handles.H2.XLim(1):2:handles.H2.XLim(2),'XTickLabel',num2str([(handles.H2.XLim(1):2:handles.H2.XLim(2))-handles.H2.XLim(1)]'))
set(handles.H2,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual');
handles.H2.XLim = tstart + timewindow;
handles.H2.FontSize = 12;
handles.H2.FontWeight = 'bold';
set(handles.H2,'YTick',[0 5],'YTickLabel',{'min' 'max'},'TickDir','out','Box','off');
handles.H2.TickLength(1) = 0.005;

% initialize another axes object to plot the timestamp bar
handles.H3 = axes;
handles.H3.Position = handles.H2.Position;
handles.H3.Position(2) = handles.H3.Position(2) - 0.005;
handles.H3.Position(4) = handles.H3.Position(4) + 0.01;
axes(handles.H3);
set(gca, 'Color', 'none');
handles.H3.YLim = handles.H2.YLim;
handles.H3.XLim = tstart + timewindow;
set(handles.H3,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.H3,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual','TickDir','out','Box','off');
handles.bar = line([tstart tstart],handles.H3.YLim,'color',[239 41 41]./256,'LineWidth',2);
handles.H3.FontSize = 14;
handles.H3.TickLength(1) = 0.005;
handles.H3.YColor = 'none';
handles.H3.XColor = 'none';

% add two more axes to display transfer function and the motor location
handles.H5 = axes;
handles.H5.Position = handles.H2.Position;
handles.H5.Position(1) = handles.H2.Position(1) + handles.H2.Position(3)+margin/20;
handles.H5.Position(3) = 0.02;
YLimfactor = handles.H2.Position(4)/(handles.H2.YLim(2) - handles.H2.YLim(1));
handles.H5.Position(2) = handles.H5.Position(2) + YLimfactor*abs(handles.H2.YLim(1));
handles.H5.Position(4) = YLimfactor*5;
handles.TF_plot = imagesc(((-50:1:50)')/50,[-1 1]);
colormap(handles.H5, brewermap([33],'rdbu'));
axis off tight
set(handles.H5,'YLim',[0 100]);

handles.H6 = axes;
handles.H6.Position = handles.H5.Position;
handles.H6.Position(1) = handles.H6.Position(1) + handles.H5.Position(3);
handles.motor_location = plot(1,2,'r<','MarkerFaceColor','k','MarkerEdgeColor','k');
axis off tight
set(handles.H6,'YLim',[0 120]);
set(handles.H6, 'Color', 'none');

% % add another axes to mark the targetzone boundaries
zonehalfwidth = 0.0625*handles.H1.Position(3);
zoneheight = 0.06;
handles.H4 = axes;
handles.H4.Position = handles.H1.Position;
handles.H4.Position(1) = handles.H1.Position(1) + (handles.H1.Position(3))/2 - zonehalfwidth;
handles.H4.Position(2) = handles.H1.Position(2) - zoneheight - 0.001;
handles.H4.Position(3) = 2*zonehalfwidth;

set(gca,'Color',[0.8 0.8 0.8]);
handles.H4.Position(3) = 2*zonehalfwidth;
handles.H4.Position(4) = zoneheight;
set(handles.H4,'XTick',[],'XTickLabel',' ','XTickMode','manual','XTickLabelMode','manual');
set(handles.H4,'YTick',[],'YTickLabel',' ','YTickMode','manual','YTickLabelMode','manual','Box','on');
handles.H4.Color = [1 1 1];

% get current TF
TrialNum = find(TrialInfo.SessionTimestamps(:,1)>=tstart,1);
MyTF = AllTFs(TrialInfo.TargetZoneType(TrialNum,:),:);
%handles.TF_plot.CData = MyTF';
handles.TF_plot.CData = flipud(MyTF')/MaxMotorLocation;

frame_offset = 0;
frameID = FrameStartIdx_Cam1 - 1;
trialstate = [0 0];
for i = 1:frames_to_show
    axes(handles.H1);
    frameID = frameID + 1;
    frame = imread(fullfile(VideoFolder{1},[filetag{1},num2str(frameID),'.png']));
    handles.frame.CData = imresize(frame,0.5);
    
    axes(handles.H1B);
    frame = imread(fullfile(VideoFolder{2},[filetag{2},num2str(frameID),'.png']));
    handles.frameB.CData = imresize(frame,0.5);
    
    %tstart = tstart + (1/framerate);
    tstart = FrameTimeStamps_Cam1(frameID);
    handles.bar.XData = [tstart tstart];
    pause(1/(framerate/slowby));
    
    % encoding of trial state
    t1 = find(abs(MyData(:,1)-tstart)==min(abs(MyData(:,1)-tstart)),1);
    if MyData(t1,6)>0
        trialstate(2) = 1;
        if mode(MyData(t1-19:t1,7))>0
            %handles.H4.Color = [0.8 0.8 0.8];
            %handles.H4.Color = [1 1 0];
            handles.H4.Color = [242 234 178]./255;
        else
            handles.H4.Color = [0.8 0.8 0.8];
            %handles.H4.Color = [239 41 41]./255;
        end
    else
        trialstate(2) = 0;
        %handles.H4.Color = [0.8 0.8 0.8];
        handles.H4.Color = [1 1 1];
    end
    
    % if trial just went off
    if (trialstate(1)-trialstate(2))==1
        TrialNum = find(TrialInfo.SessionTimestamps(:,1)>=tstart,1);
        MyTF = AllTFs(TrialInfo.TargetZoneType(TrialNum,:),:);
        handles.TF_plot.CData = MyTF';
    end
    
    % Show motorlocation
    motor_location = MyData(t1,13)/MaxMotorLocation;
    MyMap = flipud(handles.TF_plot.CData);
    [~,idx] = min(abs(MyMap-motor_location));
    handles.motor_location.YData = idx;
%     motor_location = MyData(t1,13);
%     [~,foo] = min(handles.TF_plot.CData);
%     MyMap = handles.TF_plot.CData;
%     MyMap(foo:end,1) = -1*MyMap(foo:end,1);
%     [~,idx] = min(abs(MyMap-motor_location/80));
%     idx = 100 - idx;
%     handles.motor_location.YData = idx;
    

    trialstate(1) = trialstate(2);
    
    if write_video
        myframe = getframe(h);
        writeVideo(writerObj,myframe);
    end
end

if write_video
    % cleanup
    close(writerObj);
    delete(vidObj);
    delete(writerObj);
end

end



