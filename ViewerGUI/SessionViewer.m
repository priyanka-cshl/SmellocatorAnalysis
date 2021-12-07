function varargout = SessionViewer(varargin)
% SESSIONVIEWER MATLAB code for SessionViewer.fig
%      SESSIONVIEWER, by itself, creates a new SESSIONVIEWER or raises the existing
%      singleton*.
%
%      H = SESSIONVIEWER returns the handle to a new SESSIONVIEWER or the handle to
%      the existing singleton*.
%
%      SESSIONVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SESSIONVIEWER.M with the given input arguments.
%
%      SESSIONVIEWER('Property','Value',...) creates a new SESSIONVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SessionViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SessionViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SessionViewer

% Last Modified by GUIDE v2.5 07-Dec-2021 10:10:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SessionViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @SessionViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SessionViewer is made visible.
function SessionViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SessionViewer (see VARARGIN)

% Choose default command line output for SessionViewer
handles.output = hObject;

% defaults
handles.SampleRate = 500;
handles.SessionLength.String = '100';
[Paths] = WhichComputer();
handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20210922_r0_processed.mat');
handles.TimeWindow.String = '100';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SessionViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SessionViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on slider movement.
function Scroller_Callback(hObject, eventdata, handles)
% hObject    handle to Scroller (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

newLims = handles.Scroller.Value * str2double(handles.SessionLength.String) + ...
    [0 str2double(handles.TimeWindow.String)];
set(handles.SpikesPlot,'XLim',newLims);
set(handles.PSTHPlot,'XLim',newLims);
set(handles.popPSTH,'XLim',newLims);
set(handles.BehaviorPlot,'XLim',newLims); 
set(handles.MotorPlot,'XLim',handles.SampleRate*newLims); 

% Update handles structure
guidata(hObject, handles);


function TimeWindow_Callback(hObject, eventdata, handles)
% hObject    handle to TimeWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeWindow as text
%        str2double(get(hObject,'String')) returns contents of TimeWindow as a double
Scroller_Callback(hObject, eventdata, handles);


% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.WhereSession.String)
        [WhichSession, SessionPath] = uigetfile(...
                                '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/O3/O3_20210922_r0_processed.mat',...
                                '*.mat', 'Select Behavior/Recording Session');
handles.WhereSession.String = fullfile(SessionPath,WhichSession);
end

% Load the relevant variables
load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
handles.NumUnits.String = num2str(size(SingleUnits,2));
if any(~cellfun(@isempty, TrialInfo.Perturbation))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation));
    u = unique(TrialInfo.Perturbation(x));
    handles.PerturbationList.String = u{1};
    for y = 2:size(u,1)
        handles.PerturbationList.String = [handles.PerturbationList.String,'; ',u{y}];
    end
else
    handles.PerturbationList.String = '';
end
handles.SampleRate = SampleRate;
[TracesOut] = ConcatenateTraces(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);

% create a corresponding timestamp vector
Timestamps = (1:numel(TracesOut.Lever{1}))/SampleRate;
% make the first timestamp equal to first trial minus startoffset - such
% that they match raw session timestamps
Timestamps = Timestamps + TrialInfo.SessionTimestamps(1,1) - startoffset;

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,1);
TrialStart_Ephys = TTLs.Trial(1,1);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

%% plot odor boxes on the behavior plot
axes(handles.BehaviorPlot);
hold off
for i = 1:4
    handles.(['Trial',num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.(['Trial',num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TrialInfo.SessionTimestamps((TrialInfo.Odor==i),1:2)' + TimestampAdjuster;
    if ~isempty(ValveTS)
        handles.(['Trial',num2str(i),'Plot']).Vertices = [ ...
            reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
            repmat([0 5 5 0]',size(ValveTS,2),1)];
        handles.(['Trial',num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
    end
end

% plot the target zone
handles.TargetZonePlot = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.2);
hold on;
handles.TargetZonePlot.EdgeColor = 'none';
TrialTS = TrialInfo.SessionTimestamps(:,1:2)' + TimestampAdjuster;
TZList =  TargetZones(TrialInfo.TargetZoneType,[3 1 1 3])';
handles.TargetZonePlot.Vertices = [ ...
        reshape([TrialTS(:) TrialTS(:)]', 2*numel(TrialTS), []) , ...
        TZList(:)];
handles.TargetZonePlot.Faces = reshape(1:2*numel(TrialTS),4,size(TrialTS,2))';
    
% plot the lever trace on top
plot(Timestamps + TimestampAdjuster, TracesOut.Lever{1},'k');

% plot Rewards
handles.reward_plot = plot(NaN, NaN, 'color',Plot_Colors('t'),'Linewidth',1.25);
tick_timestamps =  Timestamps(find(diff(TracesOut.Rewards{1}==1)) + 1)' + TimestampAdjuster;
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; 5; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.reward_plot,'XData',tick_x,'YData',tick_y);

set(gca,'YLim', [0 6], 'YTick', [],...
    'XTick', [], 'XLim', [0 str2double(handles.TimeWindow.String)]);

%% plot odor boxes on the spikes plot
axes(handles.SpikesPlot);
hold off
[handles] = EventsPlotter(handles,'Odor','water_plot',TTLs,TuningTTLs);

% plot all spikes
RecordingSessionOverview(SingleUnits);
set(gca,'YLim', [0 10+str2double(handles.NumUnits.String)+1], 'YTick', [],...
    'XTick', [],...
    'TickDir','out','XLim', [0 str2double(handles.TimeWindow.String)]);


%% plot odor boxes on the PSTH plot
axes(handles.PSTHPlot);
hold off
[handles] = EventsPlotter(handles,'Smell','h2o_plot',TTLs,TuningTTLs);
% plot all PSTHs
[popFR] = RecordingSessionOverview(SingleUnits,'rastermode',0,'sessionlength',str2double(handles.SessionLength.String));

set(gca,'YLim', [0 10+str2double(handles.NumUnits.String)+1], 'YTick', [],...
    'XTick', [],...
    'TickDir','out','XLim', [0 str2double(handles.TimeWindow.String)]);

%% population psth
axes(handles.popPSTH);
hold off;
[handles] = EventsPlotter(handles,'Stink','WaterPlot',TTLs,TuningTTLs);

% plot population psth
taxis = (1:size(popFR,1))/100;
plot(taxis,popFR(:,1),'color',Plot_Colors('k'));

set(gca, 'YLim', [0 ceil(max(popFR(:,1)))], ...
    'TickDir','out','XLim', [0 str2double(handles.TimeWindow.String)]);


%% update the motor plot
axes(handles.MotorPlot);
hold off
MotorTrajectory = NaN + zeros(str2double(handles.SessionLength.String)*SampleRate,1);
% Fill in the Motor Data from the behavior session period
Indices = round((Timestamps + TimestampAdjuster)*SampleRate);
MotorTrajectory(Indices,1) = TracesOut.Motor{1}/100;

% add the tuning period motor locations as well
for x = 1:size(TuningTTLs,1)
    if TuningTTLs(x,4) && ~isnan(TuningTTLs(x,5)) && ~isnan(TuningTTLs(x,7))
        Indices = round(TuningTTLs(x,4)*SampleRate):round(TuningTTLs(x,6)*SampleRate);
        MotorTrajectory(Indices,1) = TuningTTLs(x,7)/100;
    end
end
% MotorTrajectory(:,2) = MotorTrajectory(:,1);
% % keep only positive values in col 1
% MotorTrajectory((MotorTrajectory(:,1)<0),1) = 1.1;
% MotorTrajectory((MotorTrajectory(:,2)>0),2) = 1.1;
% MotorTrajectory = abs(MotorTrajectory);

alphamask = ~isnan(MotorTrajectory);

handles.MotorTrajectoryPlot = imagesc(MotorTrajectory',[-1 1]);
colormap(brewermap(100,'RdYlBu'));
set(handles.MotorTrajectoryPlot, 'AlphaData', alphamask');

set(gca, 'YTick', [], 'XTick', [], ...
    'TickDir','out','XLim', [0 SampleRate*str2double(handles.TimeWindow.String)]);

% display PSTH or Rasters as needed
PSTHView_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in NewSession.
function NewSession_Callback(hObject, eventdata, handles)
% hObject    handle to NewSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.WhereSession.String = '';


% --- Executes on button press in PSTHView.
function PSTHView_Callback(hObject, eventdata, handles)
% hObject    handle to PSTHView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PSTHView
if handles.PSTHView.Value
    set([handles.SpikesPlot; handles.SpikesPlot.Children], 'Visible','off');
    set([handles.PSTHPlot; handles.PSTHPlot.Children], 'Visible','on');
else
    set([handles.SpikesPlot; handles.SpikesPlot.Children], 'Visible','on');
    set([handles.PSTHPlot; handles.PSTHPlot.Children], 'Visible','off');
end
