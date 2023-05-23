function varargout = Smellocator_DataViewer(varargin)
% SMELLOCATOR_DATAVIEWER MATLAB code for Smellocator_DataViewer.fig
%      SMELLOCATOR_DATAVIEWER, by itself, creates a new SMELLOCATOR_DATAVIEWER or raises the existing
%      singleton*.
%
%      H = SMELLOCATOR_DATAVIEWER returns the handle to a new SMELLOCATOR_DATAVIEWER or the handle to
%      the existing singleton*.
%
%      SMELLOCATOR_DATAVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMELLOCATOR_DATAVIEWER.M with the given input arguments.
%
%      SMELLOCATOR_DATAVIEWER('Property','Value',...) creates a new SMELLOCATOR_DATAVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Smellocator_DataViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Smellocator_DataViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Smellocator_DataViewer

% Last Modified by GUIDE v2.5 22-May-2023 15:45:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Smellocator_DataViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @Smellocator_DataViewer_OutputFcn, ...
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


% --- Executes just before Smellocator_DataViewer is made visible.
function Smellocator_DataViewer_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% defaults
handles.SampleRate = 500;
handles.SessionLength.String = '100';
handles.TimeWindow.String = '20';
handles.RespirationScaling.Data = [6 0.5];
handles.WhereSession.String = '';
handles.Scroller.Value = 0.2;
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Smellocator_DataViewer_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)

global TargetZones;
global SampleRate; 
global startoffset;

if isempty(handles.WhereSession.String)
    disp('Please specify the behavior session to plot');
    return;
end

% Load the relevant variables
load(handles.WhereSession.String, 'Traces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'SampleRate');

handles.SampleRate = SampleRate;
[TracesOut, whichTraces] = ConcatenateTraces2Matrix(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);
% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(whichTraces,'Timestamps')));

%% Calculate Trial On-Off timestamps
TrialColumn = TracesOut(:,find(strcmp(whichTraces,'TrialState')));
TrialColumn(TrialColumn~=0) = 1; % make logical
TrialOn = find(diff([0; TrialColumn])>0);
TrialOff =  find(diff(TrialColumn)<0)+1;

% account for cases where acquisition started/ended in between a trial
while TrialOn(1)>TrialOff(1)
    TrialOff(1,:) = [];
end
while TrialOn(end)>TrialOff(end)
    TrialOn(end,:) = [];
end 

TrialIndices = [TrialOn TrialOff];
TrialTimeStamps = [Timestamps(TrialOn,1) Timestamps(TrialOff,1)];

handles.SessionLength.String = TrialTimeStamps(end,2);

%% plot odor boxes on the behavior plot
axes(handles.BehaviorPlot);
hold off
for i = 1:4
    handles.(['Trial',num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.(['Trial',num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TrialTimeStamps((TrialInfo.Odor==i),1:2)';
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
TrialTS = TrialTimeStamps(:,1:2)';
TZList =  TargetZones(TrialInfo.TargetZoneType,[3 1 1 3])';
handles.TargetZonePlot.Vertices = [ ...
        reshape([TrialTS(:) TrialTS(:)]', 2*numel(TrialTS), []) , ...
        TZList(:)];
handles.TargetZonePlot.Faces = reshape(1:2*numel(TrialTS),4,size(TrialTS,2))';
    
% plot the lever trace on top
plot(Timestamps, TracesOut(:,find(strcmp(whichTraces,'Lever'))),'k');

% plot Rewards
handles.reward_plot = plot(NaN, NaN, 'color',Plot_Colors('t'),'Linewidth',1.25);
tick_timestamps =  Timestamps(find(diff(TracesOut(:,find(strcmp(whichTraces,'Rewards')))==1)) + 1);
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; 5; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.reward_plot,'XData',tick_x,'YData',tick_y);

% plot Licks
handles.lick_plot = plot(NaN, NaN, 'color',Plot_Colors('r'),'Linewidth',1.25);
tick_timestamps =  Timestamps(find(diff(TracesOut(:,find(strcmp(whichTraces,'Licks')))==1)) + 1);
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [5.2; 5.5; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.lick_plot,'XData',tick_x,'YData',tick_y);

set(gca,'YLim', [0 7], 'YTick', [],...
    'XTick', [], 'XLim', [0 str2double(handles.TimeWindow.String)]);

% plot respiration
handles.SniffData = TracesOut(:,find(strcmp(whichTraces,'Sniffs')));
handles.SniffData = handles.SniffData - mean(handles.SniffData);

handles.respiration_plot = plot(Timestamps, ...
    handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.SniffData, ...
    'color',Plot_Colors('b'),'Linewidth',1);

%% update the motor plot
axes(handles.MotorPlot);
hold off
MotorTrajectory = NaN + zeros(str2double(handles.SessionLength.String)*SampleRate,1);
% Fill in the Motor Data from the behavior session period
Indices = round((Timestamps)*SampleRate);
MotorTrajectory(Indices,1) = TracesOut(:,find(strcmp(whichTraces,'OdorLocation')))/100;

alphamask = ~isnan(MotorTrajectory);

handles.MotorTrajectoryPlot = imagesc(MotorTrajectory',[-1 1]);
colormap(brewermap(100,'RdYlBu'));
set(handles.MotorTrajectoryPlot, 'AlphaData', alphamask');

set(gca, 'YTick', [], 'XTick', [], ...
    'TickDir','out','XLim', [0 SampleRate*str2double(handles.TimeWindow.String)]);

Scroller_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in NewSession.
function NewSession_Callback(hObject, eventdata, handles)
handles.WhereSession.String = '';

% --- Executes on slider movement.
function Scroller_Callback(hObject, eventdata, handles)
newLims = handles.Scroller.Value * str2double(handles.SessionLength.String) + ...
    [0 str2double(handles.TimeWindow.String)];
set(handles.BehaviorPlot,'XLim',newLims); 
set(handles.MotorPlot,'XLim',handles.SampleRate*newLims); 

% Update handles structure
guidata(hObject, handles);

function TimeWindow_Callback(hObject, eventdata, handles)
Scroller_Callback(hObject, eventdata, handles);

% --- Executes when entered data in editable cell(s) in RespirationScaling.
function RespirationScaling_CellEditCallback(hObject, eventdata, handles)
handles.respiration_plot.YData = handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.SniffData';
guidata(hObject, handles);
