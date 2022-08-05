function varargout = SessionViewerBasic(varargin)
% SESSIONVIEWERBASIC MATLAB code for SessionViewerBasic.fig
%      SESSIONVIEWERBASIC, by itself, creates a new SESSIONVIEWERBASIC or raises the existing
%      singleton*.
%
%      H = SESSIONVIEWERBASIC returns the handle to a new SESSIONVIEWERBASIC or the handle to
%      the existing singleton*.
%
%      SESSIONVIEWERBASIC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SESSIONVIEWERBASIC.M with the given input arguments.
%
%      SESSIONVIEWERBASIC('Property','Value',...) creates a new SESSIONVIEWERBASIC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SessionViewerBasic_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SessionViewerBasic_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SessionViewerBasic

% Last Modified by GUIDE v2.5 28-Jul-2022 08:43:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SessionViewerBasic_OpeningFcn, ...
                   'gui_OutputFcn',  @SessionViewerBasic_OutputFcn, ...
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


% --- Executes just before SessionViewerBasic is made visible.
function SessionViewerBasic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SessionViewerBasic (see VARARGIN)

% Choose default command line output for SessionViewerBasic
handles.output = hObject;

% defaults
handles.SampleRate = 500;
handles.SessionLength.String = '100';
[Paths] = WhichComputer();
handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20210922_r0_processed.mat');
handles.TimeWindow.String = '20';
handles.RespirationScaling.Data = [6 0.5];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SessionViewerBasic wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SessionViewerBasic_OutputFcn(hObject, eventdata, handles) 
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
global TargetZones;

if isempty(handles.WhereSession.String)
    [Paths] = WhichComputer();
    if ~isempty(handles.MouseName.String)
        DefaultPath = [handles.MouseName.String,filesep,handles.MouseName.String,'_*.mat'];
    else
        DefaultPath = 'O9/O9_20220625_r0.mat';
    end
    [WhichSession, SessionPath] = uigetfile(...
                                fullfile(Paths.Grid.Behavior,DefaultPath),...
                                'Select Behavior Session');
    handles.WhereSession.String = fullfile(SessionPath,WhichSession);
end

% Load the relevant variables
if isempty(strfind(handles.WhereSession.String,'_processed.mat'))
    [MyData, MySettings, DataTags] = ReadSessionData(handles.WhereSession.String);
    [Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
    [MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags);
    [Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);
    TTLs = []; SampleRate = 500; startoffset = 1;
else    
    load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');
end

if ~isempty(TTLs)
    handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
    handles.NumUnits.String = num2str(size(SingleUnits,2));
else
    handles.SessionLength.String = TrialInfo.SessionTimestamps(end,2);
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
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);

if ~isempty(TTLs)
    TrialStart_Ephys = TTLs.Trial(1,2);
    % factor to convert all behavior timestamps to match Ephys
    TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;
else
    TrialStart_Ephys = 0;
    TimestampAdjuster = 0;
end

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
%     if any(strcmp(TrialInfo.Perturbation(x),'Halt-Flip')) || ...
%             any(strcmp(TrialInfo.Perturbation(x),'NoOdor')) || ...
%                 any(strcmp(TrialInfo.Perturbation(x),'Offset-II')) || ...
%                     any(strcmp(TrialInfo.Perturbation(x),'Halt-II'))
%         PerturbTS = TrialInfo.SessionTimestamps(x,1:2)' + TimestampAdjuster;
%     end
    
    if any(strcmp(TrialInfo.Perturbation(x),'RuleReversal'))
        y = find(strcmp(TrialInfo.Perturbation(:,1),'RuleReversal'));
        PerturbTS = TrialInfo.SessionTimestamps(y,1:2)' + TimestampAdjuster;
    else % 'Halt-Flip' 'NoOdor' 'Offset-II' 'Halt-II' 'GainChange'
        PerturbTS = TrialInfo.SessionTimestamps(x,1:2)' + TimestampAdjuster;
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

set(gca,'YLim', [0 7], 'YTick', [],...
    'XTick', [], 'XLim', [0 str2double(handles.TimeWindow.String)]);


% plot respiration
handles.SniffData = TracesOut.Sniffs{1};
handles.SniffData = handles.SniffData - mean(handles.SniffData);

handles.respiration_plot = plot(Timestamps + TimestampAdjuster, ...
    handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.SniffData, ...
    'color',Plot_Colors('r'),'Linewidth',1);

%% update the motor plot
axes(handles.MotorPlot);
hold off
MotorTrajectory = NaN + zeros(str2double(handles.SessionLength.String)*SampleRate,1);
% Fill in the Motor Data from the behavior session period
Indices = round((Timestamps + TimestampAdjuster)*SampleRate);
MotorTrajectory(Indices,1) = TracesOut.Motor{1}/100;

alphamask = ~isnan(MotorTrajectory);

handles.MotorTrajectoryPlot = imagesc(MotorTrajectory',[-1 1]);
colormap(brewermap(100,'RdYlBu'));
set(handles.MotorTrajectoryPlot, 'AlphaData', alphamask');

set(gca, 'YTick', [], 'XTick', [], ...
    'TickDir','out','XLim', [0 SampleRate*str2double(handles.TimeWindow.String)]);

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


% --- Executes on button press in ScrollUp.
function ScrollUp_Callback(hObject, eventdata, handles)
% hObject    handle to ScrollUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currLim = handles.SpikesPlot.YLim;
newLim = currLim + round(diff(currLim)/2);
% newLim(1) = max(newLim(1), 0);
% newLim(2) = min(newLim(2),(10+str2double(handles.NumUnits.String)+1));
set(handles.SpikesPlot,'YLim', newLim);
set(handles.PSTHPlot,'YLim', newLim);


% --- Executes on button press in ScrollDown.
function ScrollDown_Callback(hObject, eventdata, handles)
% hObject    handle to ScrollDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currLim = handles.SpikesPlot.YLim;
newLim = currLim - round(diff(currLim)/2);
% newLim(1) = max(newLim(1), 0);
% newLim(2) = min(newLim(2),(10+str2double(handles.NumUnits.String)+1));
set(handles.SpikesPlot,'YLim', newLim);
set(handles.PSTHPlot,'YLim', newLim);


function ZoomVal_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZoomVal as text
%        str2double(get(hObject,'String')) returns contents of ZoomVal as a double
zoomfactor = str2double(handles.ZoomVal.String);
OriginalLim = 10+str2double(handles.NumUnits.String)+1;
newLim = [0 ceil(OriginalLim/zoomfactor)];
set(handles.SpikesPlot,'YLim', newLim);
set(handles.PSTHPlot,'YLim', newLim);



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


% --- Executes when entered data in editable cell(s) in RespirationScaling.
function RespirationScaling_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to RespirationScaling (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.respiration_plot.YData = handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.SniffData';

% Update handles structure
guidata(hObject, handles);



function MouseName_Callback(hObject, eventdata, handles)
% hObject    handle to MouseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MouseName as text
%        str2double(get(hObject,'String')) returns contents of MouseName as a double


% --- Executes during object creation, after setting all properties.
function MouseName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MouseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
