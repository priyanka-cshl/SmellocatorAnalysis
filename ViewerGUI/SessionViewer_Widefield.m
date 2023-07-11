function varargout = SessionViewer_Widefield(varargin)
% SESSIONVIEWER_WIDEFIELD MATLAB code for SessionViewer_Widefield.fig
%      SESSIONVIEWER_WIDEFIELD, by itself, creates a new SESSIONVIEWER_WIDEFIELD or raises the existing
%      singleton*.
%
%      H = SESSIONVIEWER_WIDEFIELD returns the handle to a new SESSIONVIEWER_WIDEFIELD or the handle to
%      the existing singleton*.
%
%      SESSIONVIEWER_WIDEFIELD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SESSIONVIEWER_WIDEFIELD.M with the given input arguments.
%
%      SESSIONVIEWER_WIDEFIELD('Property','Value',...) creates a new SESSIONVIEWER_WIDEFIELD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SessionViewer_Widefield_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SessionViewer_Widefield_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SessionViewer_Widefield

% Last Modified by GUIDE v2.5 11-Jul-2023 07:48:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SessionViewer_Widefield_OpeningFcn, ...
                   'gui_OutputFcn',  @SessionViewer_Widefield_OutputFcn, ...
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


% --- Executes just before SessionViewer_Widefield is made visible.
function SessionViewer_Widefield_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SessionViewer_Widefield (see VARARGIN)

% Choose default command line output for SessionViewer_Widefield
handles.output = hObject;

% defaults
handles.SampleRate = 500;
handles.SessionLength.String = '100';
[Paths] = WhichComputer();
handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20210922_r0_processed.mat');
handles.TimeWindow.String = '100';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SessionViewer_Widefield wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SessionViewer_Widefield_OutputFcn(hObject, eventdata, handles) 
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
    [Paths] = WhichComputer();
        [WhichSession, SessionPath] = uigetfile(...
                                fullfile(Paths.ProcessedSessions,'HX3/HX3_20230505_r0_processed.mat'),...
                                'Select Behavior or Recording Session');
handles.WhereSession.String = fullfile(SessionPath,WhichSession);
end

% Load the relevant variables
load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');
           
load(handles.ROIPath.String, 'WhichROIs','refPix','C','stackdims');           

handles.SessionLength.String = TrialInfo.SessionTimestamps(end,2);
if ~isempty(WhichROIs)
    handles.NumUnits.String = num2str(size(WhichROIs,1));
else
    handles.NumUnits.String = 'NaN';
end

if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    u = unique(TrialInfo.Perturbation(x));
    handles.PerturbationList.String = u{1};
    for y = 2:size(u,1)
        handles.PerturbationList.String = [handles.PerturbationList.String,'; ',u{y}];
    end
else
    x = [];
    handles.PerturbationList.String = '';
end
handles.SampleRate = SampleRate;
[TracesOut] = ConcatenateTraces(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);

% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% calculate the timestamp difference between Ephys and Behavior
TimestampAdjuster = 0; % by definition

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

% plot the perturbation periods
handles.(['Perturbation',num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors('r'));
hold on;
handles.(['Perturbation',num2str(i),'Plot']).EdgeColor = 'none';


if ~isempty(x)
    % for replay sessions
    if any(strcmp(TrialInfo.Perturbation(x),'OL-Template'))
        % get start and stop TS of the template
        templateStart = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'first'));
        templateEnd   = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'last'));
        PerturbTS(1,1) = TrialInfo.SessionTimestamps(templateStart,1) + TimestampAdjuster;
        PerturbTS(2,1) = TrialInfo.SessionTimestamps(templateEnd,1) + TimestampAdjuster;
        % get start and stop of the replays
        replays = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Replay')));
        PerturbTS = horzcat(PerturbTS, ...
            TrialInfo.SessionTimestamps(replays,1:2)' + TimestampAdjuster );
    end
    
    % for other perturbations
    if any(strcmp(TrialInfo.Perturbation(x),'Halt-Flip'))
        PerturbTS = TrialInfo.SessionTimestamps(x,1:2)' + TimestampAdjuster;
    end
    
    if any(strcmp(TrialInfo.Perturbation(x),'RuleReversal'))
        y = find(strcmp(TrialInfo.Perturbation(:,1),'RuleReversal'));
        PerturbTS = TrialInfo.SessionTimestamps(y,1:2)' + TimestampAdjuster;
    end
    
    if ~isempty(PerturbTS)
        handles.(['Perturbation',num2str(i),'Plot']).Vertices = [ ...
            reshape([PerturbTS(:) PerturbTS(:)]', 2*numel(PerturbTS), []) , ...
            repmat([5 5.5 5.5 5]',size(PerturbTS,2),1)];
        handles.(['Perturbation',num2str(i),'Plot']).Faces = reshape(1:2*numel(PerturbTS),4,size(PerturbTS,2))';
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

if ~isempty(WhichROIs)    
    
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
    
end

%% update the motor plot
axes(handles.MotorPlot);
hold off
MotorTrajectory = NaN + zeros(str2double(handles.SessionLength.String)*SampleRate,1);
% Fill in the Motor Data from the behavior session period
Indices = round((Timestamps + TimestampAdjuster)*SampleRate);
MotorTrajectory(Indices,1) = TracesOut.Motor{1}/100;

if ~isempty(TuningTTLs)
    
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
end
alphamask = ~isnan(MotorTrajectory);

handles.MotorTrajectoryPlot = imagesc(MotorTrajectory',[-1 1]);
colormap(brewermap(100,'RdYlBu'));
set(handles.MotorTrajectoryPlot, 'AlphaData', alphamask');

set(gca, 'YTick', [], 'XTick', [], ...
    'TickDir','out','XLim', [0 SampleRate*str2double(handles.TimeWindow.String)]);

if ~isempty(TTLs)
    % display PSTH or Rasters as needed
    PSTHView_Callback(hObject, eventdata, handles);
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in NewSession.
function NewSession_Callback(hObject, eventdata, handles)
% hObject    handle to NewSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.WhereSession.String = '';
handles.ROIPath.String = '';


% --- Executes on button press in PSTHView.

function WhichUnits_Callback(hObject, eventdata, handles)
% hObject    handle to WhichUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhichUnits as text
%        str2double(get(hObject,'String')) returns contents of WhichUnits as a double


% --- Executes during object creation, after setting all properties.
function WhichUnits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhichUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
