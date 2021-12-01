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

% Last Modified by GUIDE v2.5 01-Dec-2021 16:55:38

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
handles.SessionLength = 100;

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

newLims = get(hObject,'Value')*handles.SessionLength + [0 str2double(handles.TimeWindow.String)];
set(handles.SpikesPlot,'XLim',newLims);
set(handles.BehaviorPlot,'XLim',newLims); 

% --- Executes during object creation, after setting all properties.
function Scroller_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scroller (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function TimeWindow_Callback(hObject, eventdata, handles)
% hObject    handle to TimeWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeWindow as text
%        str2double(get(hObject,'String')) returns contents of TimeWindow as a double
Scroller_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function TimeWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[WhichSession, SessionPath] = uigetfile(...
                                '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/O5/O5_20210923_r0_processed.mat',...
                                '*.mat', 'Select Behavior/Recording Session');
                            
% Load the relevant variables
load(fullfile(SessionPath,WhichSession), 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength = 10*ceil(TTLs.Trial(end,2)/10);
handles.TotalUnits = size(SingleUnits,2);
[TracesOut] = ConcatenateTraces(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);
TrialStart_behavior = TrialInfo.SessionTimestamps(1,1);
TrialStart_Ephys = TTLs.Trial(1,1);
% convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

% plot odor boxes on the behavior plot
axes(handles.BehaviorPlot);
for i = 1:3
    handles.(['Trial',num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.(['Trial',num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TrialInfo.SessionTimestamps((TrialInfo.Odor==i),1:2)' + TimestampAdjuster;
    handles.(['Trial',num2str(i),'Plot']).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat([0 5 5 0]',size(ValveTS,2),1)];
    handles.(['Trial',num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end
% plot the lever trace on top
Timestamps = (1:numel(TracesOut.Lever{1}))/SampleRate;
% make the first timestamp equal to first trial minus startoffset
Timestamps = Timestamps - Timestamps(1) + TrialInfo.SessionTimestamps(1,1) - startoffset;
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

set(gca,'YLim', [0 8], 'YTick', [],...
    'XTick', [], 'XLim', [0 str2double(handles.TimeWindow.String)]);

% plot odor boxes on the spikes plot
axes(handles.SpikesPlot);
for i = 1:3
    handles.(['Odor',num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.(['Odor',num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TTLs.(['Odor',num2str(i)])(:,1:2)';
    handles.(['Odor',num2str(i),'Plot']).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat(handles.TotalUnits*[0 1 1 0]',size(ValveTS,2),1)];
    handles.(['Odor',num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end

% plot Rewards
handles.water_plot = plot(NaN, NaN, 'color',Plot_Colors('t'),'Linewidth',1.25);
tick_timestamps =  TTLs.Reward(:,1);
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; handles.TotalUnits; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.water_plot,'XData',tick_x,'YData',tick_y);

% plot all spikes
RecordingSessionOverview(SingleUnits);
set(gca,'YLim', [0 handles.TotalUnits], 'YTick', [],...
    'TickDir','out','XLim', [0 str2double(handles.TimeWindow.String)]);
% Update handles structure
guidata(hObject, handles);
