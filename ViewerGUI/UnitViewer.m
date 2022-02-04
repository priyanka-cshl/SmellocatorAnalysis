function varargout = UnitViewer(varargin)
% UNITVIEWER MATLAB code for UnitViewer.fig
%      UNITVIEWER, by itself, creates a new UNITVIEWER or raises the existing
%      singleton*.
%
%      H = UNITVIEWER returns the handle to a new UNITVIEWER or the handle to
%      the existing singleton*.
%
%      UNITVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNITVIEWER.M with the given input arguments.
%
%      UNITVIEWER('Property','Value',...) creates a new UNITVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UnitViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UnitViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UnitViewer

% Last Modified by GUIDE v2.5 04-Feb-2022 16:50:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UnitViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @UnitViewer_OutputFcn, ...
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


% --- Executes just before UnitViewer is made visible.
function UnitViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UnitViewer (see VARARGIN)

% Choose default command line output for UnitViewer
handles.output = hObject;

[Paths] = WhichComputer();
handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
handles.CurrentUnit.Data(1) = NaN;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UnitViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UnitViewer_OutputFcn(hObject, eventdata, handles) 
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
if isempty(handles.WhereSession.String)
    [Paths] = WhichComputer();
        [WhichSession, SessionPath] = uigetfile(...
                                fullfile(Paths.ProcessedSessions,'O3/O3_20210922_r0_processed.mat'),...
                                'Select Behavior or Recording Session');
handles.WhereSession.String = fullfile(SessionPath,WhichSession);
end

%% get the data loaded
MySession = handles.WhereSession.String;
[TracesOut, ColNames, handles.TrialInfo, SingleUnits, TTLs, SampleRate, TimestampAdjuster] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

%% Get all spikes, all units aligned to trials
[handles.AlignedSpikes, handles.Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(handles.TrialInfo.TrialID,2),handles.TrialInfo);

handles.NumUnits.String = num2str(size(SingleUnits,2));
handles.CurrentUnit.Data(1) = 1;

UpdatePlots(handles);

% Update handles structure
guidata(hObject, handles);

function UpdatePlots(handles)
whichUnit = handles.CurrentUnit.Data(1);
AlignType = handles.AlignTo.Value;
for i = 1:3
    axes(handles.(['axes',num2str(i)])); 
    cla reset; 
    hold on
    PlotFullSession(whichUnit, i, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignType);
    switch AlignType
        case {1,2}
            set(gca, 'XLim', [-1.2 6]);
        case {3,4}
            set(gca, 'XLim', [-5.2 1]);
        case 5
            set(gca, 'XLim', [-2.2 5]);
    end
end

function WhereSession_Callback(hObject, eventdata, handles)
% hObject    handle to WhereSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhereSession as text
%        str2double(get(hObject,'String')) returns contents of WhereSession as a double


% --- Executes on button press in NewSession.
function NewSession_Callback(hObject, eventdata, handles)
% hObject    handle to NewSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in NextUnit.
function NextUnit_Callback(hObject, eventdata, handles)
% hObject    handle to NextUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CurrentUnit.Data(1) = min(handles.CurrentUnit.Data(1)+1,str2double(handles.NumUnits.String));
UpdatePlots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in PrevUnit.
function PrevUnit_Callback(hObject, eventdata, handles)
% hObject    handle to PrevUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CurrentUnit.Data(1) = max(handles.CurrentUnit.Data(1)-1,1);
UpdatePlots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in AlignTo.
function AlignTo_Callback(hObject, eventdata, handles)
% hObject    handle to AlignTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdatePlots(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns AlignTo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AlignTo


% --- Executes when entered data in editable cell(s) in CurrentUnit.
function CurrentUnit_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to CurrentUnit (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
UpdatePlots(handles);
% Update handles structure
guidata(hObject, handles);
