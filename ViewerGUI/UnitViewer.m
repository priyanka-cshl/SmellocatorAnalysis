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

% Last Modified by GUIDE v2.5 11-Mar-2022 13:05:53

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
[TracesOut, ColNames, handles.TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

%% Get all spikes, all units aligned to trials
[handles.AlignedSpikes, handles.Events, handles.whichtetrode] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession);

if any(strcmp(handles.TrialInfo.Perturbation,'OL-Replay'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events);
end
handles.NumUnits.String = num2str(size(SingleUnits,2));
handles.CurrentUnit.Data(1) = 1;

UpdatePlots(handles);

% Update handles structure
guidata(hObject, handles);

function UpdatePlots(handles)
whichUnit = handles.CurrentUnit.Data(1);
AlignType = handles.AlignTo.Value;
MyColors1 = brewermap(15,'*PuBu');
MyColors2 = brewermap(15,'*OrRd');
handles.tetrode.String = num2str(handles.whichtetrode(whichUnit));
switch AlignType
    case {1,2,6}
        myXlim = [-1.2 6];
    case {3,4}
        myXlim = [-5.2 1];
    case 5
        myXlim = [-2.2 5];
end
    
for i = 1:3
    axes(handles.(['axes',num2str(i)])); 
    cla reset; 
    hold on
    % plot baseline trials
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'RuleReversal'))
        [trialsdone, FRs, BinOffset, P_FRs] = PlotRuleReversalSession(whichUnit, i, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignType);
    else
        [trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, i, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignType);
    end
    if any(strcmp(handles.TrialInfo.Perturbation,'OL-Replay'))
        % plot replay trials
        AddReplay2FullSession(trialsdone, whichUnit, i, handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo, AlignType, handles.SortReplay.Value);
    end
    set(gca, 'XLim', myXlim);
    
    axes(handles.(['axes',num2str(i+3)])); 
    cla reset; 
    hold on
    
    if ~any(strcmp(handles.TrialInfo.Perturbation(:,1),'RuleReversal'))
        for t = 1:size(FRs,1)
            plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(t,:),'Color',MyColors1(t,:),'Linewidth',1);
        end
        
        if ~handles.HidePSTH2.Value
            if ~isempty(P_FRs)
                for t = 1:size(FRs,1)
                    set(groot,'defaultAxesColorOrder',MyColors2);
                    plot((1:size(P_FRs,2))*0.002+BinOffset/1000,P_FRs(t,:),'Color',MyColors2(t,:),'Linewidth',1);
                end
            end
        end
    end
    set(gca, 'XLim', myXlim);
    
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
handles.WhereSession.String = [];
handles.NumUnits.String = '';
handles.CurrentUnit.Data(1) = NaN;
% Update handles structure
guidata(hObject, handles);


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


% --- Executes on button press in SortReplay.
function SortReplay_Callback(hObject, eventdata, handles)
% hObject    handle to SortReplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SortReplay
UpdatePlots(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in HidePSTH2.
function HidePSTH2_Callback(hObject, eventdata, handles)
% hObject    handle to HidePSTH2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdatePlots(handles);
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of HidePSTH2



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to tetrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tetrode as text
%        str2double(get(hObject,'String')) returns contents of tetrode as a double


% --- Executes during object creation, after setting all properties.
function tetrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tetrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
