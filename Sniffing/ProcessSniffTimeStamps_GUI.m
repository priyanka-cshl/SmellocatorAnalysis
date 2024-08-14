function varargout = ProcessSniffTimeStamps_GUI(varargin)
% PROCESSSNIFFTIMESTAMPS_GUI MATLAB code for ProcessSniffTimeStamps_GUI.fig
%      PROCESSSNIFFTIMESTAMPS_GUI, by itself, creates a new PROCESSSNIFFTIMESTAMPS_GUI or raises the existing
%      singleton*.
%
%      H = PROCESSSNIFFTIMESTAMPS_GUI returns the handle to a new PROCESSSNIFFTIMESTAMPS_GUI or the handle to
%      the existing singleton*.
%
%      PROCESSSNIFFTIMESTAMPS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCESSSNIFFTIMESTAMPS_GUI.M with the given input arguments.
%
%      PROCESSSNIFFTIMESTAMPS_GUI('Property','Value',...) creates a new PROCESSSNIFFTIMESTAMPS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProcessSniffTimeStamps_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProcessSniffTimeStamps_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProcessSniffTimeStamps_GUI

% Last Modified by GUIDE v2.5 14-Aug-2024 17:31:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcessSniffTimeStamps_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcessSniffTimeStamps_GUI_OutputFcn, ...
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


% --- Executes just before ProcessSniffTimeStamps_GUI is made visible.
function ProcessSniffTimeStamps_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProcessSniffTimeStamps_GUI (see VARARGIN)

% Choose default command line output for ProcessSniffTimeStamps_GUI
handles.output = hObject;

% some initializations
handles.SniffTrace.Raw          = [];
handles.SniffTrace.Filtered     = [];
handles.SniffTrace.Timestamps   = [];
handles.OdorLocationTrace       = [];
handles.SniffsTS                = [];
handles.SessionLength           = [];

% For loading the processed session
[Paths] = WhichComputer();

if ~isempty(varargin)
    if exist(varargin{1}) == 2
        handles.WhereSession.String = varargin{1};
    else
        MouseName = regexprep(varargin{1},'_(\w+)_processed.mat','');
        handles.WhereSession.String = fullfile(Paths.ProcessedSessions,MouseName,varargin{1});
    end
else
    handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
end

% Update handles structure
guidata(hObject, handles);

if exist(handles.WhereSession.String)==2
    LoadSession_Callback(hObject, eventdata, handles);
end

% UIWAIT makes ProcessSniffTimeStamps_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ProcessSniffTimeStamps_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check that a valid session is specified
if isempty(handles.WhereSession.String)
    [Paths] = WhichComputer();
        [WhichSession, SessionPath] = uigetfile(...
                                fullfile(Paths.ProcessedSessions,'O3/O3_20210922_r0_processed.mat'),...
                                'Select Behavior or Recording Session');
handles.WhereSession.String = fullfile(SessionPath,WhichSession);
end

% load the relevant data
load(handles.WhereSession.String,'Traces','TrialInfo','SampleRate');
handles.SessionLength = ceil(Traces.Timestamps{end}(end));

% reprocess sniff traces on a trialwise basis
Traces.OdorLocation     = Traces.Motor;
[SniffTimeStamps] = ...
    TrialWiseSniffs(TrialInfo,Traces); % [sniffstart sniffstop nextsniff odorlocation sniffslope stimstate trialID]
% remove overlapping sniffs
handles.SniffsTS = SniffTimeStamps(find(SniffTimeStamps(:,end)>0),:);

% Make long concatenated traces for plotting
[TracesOut, whichtraces] = ConcatenateTraces2Mat(Traces);
handles.SniffTrace.Timestamps   = TracesOut(:,find(strcmp(whichtraces,'Timestamps')));
handles.OdorLocationTrace       = TracesOut(:,find(strcmp(whichtraces,'Motor')));
handles.SniffTrace.Raw          = TracesOut(:,find(strcmp(whichtraces,'Sniffs')));
handles.SniffTrace.Filtered     = FilterThermistor(handles.SniffTrace.Raw);

% find trace indices that correspond to detected sniff timestamps
for n = 1:size(handles.SniffsTS,1)
    % inhalation start
    [~,idx] = min(abs(handles.SniffTrace.Timestamps - handles.SniffsTS(n,1)));
    if abs(handles.SniffTrace.Timestamps(idx) - handles.SniffsTS(n,1)) < 0.004
        handles.SniffsTS(n,8) = idx;
    end
    % inhalation end
    [~,idx] = min(abs(handles.SniffTrace.Timestamps - handles.SniffsTS(n,2)));
    if abs(handles.SniffTrace.Timestamps(idx) - handles.SniffsTS(n,2)) < 0.004
        handles.SniffsTS(n,9) = idx;
    end
end

% plot both filtered and raw traces, and detcted timestamps
%set(handles.SniffingFiltered, 'nextplot', 'add');
axes(handles.SniffingFiltered);
set(handles.SniffingFiltered, 'ButtonDownFcn', {@Click2Callback, handles, hObject})
hold all
plot(handles.SniffTrace.Timestamps,handles.SniffTrace.Filtered,'color',Plot_Colors('b'));
zoom off
%hold on
% overlay detected timestamps on the plot
plot(handles.SniffTrace.Timestamps(handles.SniffsTS(:,8)),handles.SniffTrace.Filtered(handles.SniffsTS(:,8)),'ok');

plot(handles.SniffTrace.Timestamps(handles.SniffsTS(:,9)),handles.SniffTrace.Filtered(handles.SniffsTS(:,9)),'or'); %,% ...
    %'ButtonDownFcn', {@Click2Callback, handles, hObject});
set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)], ...
    'YTick', []);
currlims = get(gca,'YLim');
% Update handles structure
guidata(hObject, handles);

% Trial times
Xvals = [TrialInfo.SessionTimestamps(:,1) TrialInfo.SessionTimestamps(:,1) nan(size(TrialInfo.SessionTimestamps,1),1)]';
YVals = repmat([2*currlims NaN],size(TrialInfo.SessionTimestamps,1),1)';
plot(Xvals(:),YVals(:),'--','color',Plot_Colors('pd'));

Xvals = [TrialInfo.SessionTimestamps(:,2) TrialInfo.SessionTimestamps(:,2) nan(size(TrialInfo.SessionTimestamps,1),1)]';
YVals = repmat([2*currlims NaN],size(TrialInfo.SessionTimestamps,1),1)';
plot(Xvals(:),YVals(:),'-.','color',Plot_Colors('pl'));

set(gca,'YLim',currlims);

axes(handles.SniffingRaw);
zoom off
hold off
plot(handles.SniffTrace.Timestamps,handles.SniffTrace.Raw,'color',Plot_Colors('b'));
hold on
% overlay detected timestamps on the plot
plot(handles.SniffTrace.Timestamps(handles.SniffsTS(:,8)),handles.SniffTrace.Raw(handles.SniffsTS(:,8)),'ok');
plot(handles.SniffTrace.Timestamps(handles.SniffsTS(:,9)),handles.SniffTrace.Raw(handles.SniffsTS(:,9)),'or');
set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)], ...
     'YTick', [],  'XTick', []);


axes(handles.SniffingFiltered);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ZoomON.
function ZoomON_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom yon
guidata(hObject, handles);

% --- Executes on button press in ZoomOFF.
function ZoomOFF_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomOFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom off
guidata(hObject, handles);

% --- Executes on button press in Flag_Stretch.
function Flag_Stretch_Callback(hObject, eventdata, handles)
% hObject    handle to Flag_Stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on slider movement.
function Scroller_Callback(hObject, eventdata, handles)
% hObject    handle to Scroller (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newLims = handles.Scroller.Value * handles.SessionLength + ...
    [0 str2double(handles.WindowSize.String)];
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);

% Update handles structure
guidata(hObject, handles);


function WindowSize_Callback(hObject, eventdata, handles)
% hObject    handle to WindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowSize as text
%        str2double(get(hObject,'String')) returns contents of WindowSize as a double




% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function SniffingFiltered_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to SniffingFiltered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function SniffingFiltered_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to SniffingFiltered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coordinates = round(get(gca,'CurrentPoint'));


% --- Executes on button press in AddValleys.
function AddValleys_Callback(hObject, eventdata, handles)
% hObject    handle to AddValleys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom off
[x,y] = ginput;
% for every chosen valley
newTS = [];
for n = 1:numel(x)
    % find the trace idx
    [~,idx] = min(abs(handles.SniffTrace.Timestamps - x(n)));
    % take a window 30 ms on either side
    [valleyval,valleyidx] = min(handles.SniffTrace.Filtered(idx+[-15:1:15]));
    plot(handles.SniffTrace.Timestamps(idx+valleyidx-15),valleyval,'og');
    handles.newValleys(n,:) = [handles.SniffTrace.Timestamps(idx+valleyidx-15) valleyval];
end
uiwait(handles.figure1);
guidata(hObject, handles);

% --- Executes on button press in AddPeaks.
function AddPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to AddPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom off
[x,y] = ginput;
% for every chosen peak
for n = 1:numel(x)
    % find the trace idx
    [~,idx] = min(abs(handles.SniffTrace.Timestamps - x(n)));
    % take a window 30 ms on either side
    [valleyval,valleyidx] = max(handles.SniffTrace.Filtered(idx+[-15:1:15]));
    plot(valleyidx,valleyval,'og');
end
keyboard;
guidata(hObject, handles);

% --- Executes on button press in RemovePoints.
function RemovePoints_Callback(hObject, eventdata, handles)
% hObject    handle to RemovePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom off
[x,y] = ginput;
% for every chosen valley
for n = 1:numel(x)
    [~,idx] = min(abs(handles.SniffsTS(:,1) - x(n)));



    % find the trace idx
    [~,idx] = min(abs(handles.SniffTrace.Timestamps - x(n)));
    % take a window 30 ms on either side
    [valleyval,valleyidx] = min(handles.SniffTrace.Filtered(idx+[-15:1:15]));
    plot(valleyidx,valleyval,'og');
end
keyboard;


% --- Executes on button press in RedoStretch.
function RedoStretch_Callback(hObject, eventdata, handles)
% hObject    handle to RedoStretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi = drawrectangle;
keyboard;
