function varargout = TuningViewer(varargin)
% TUNINGVIEWER MATLAB code for TuningViewer.fig
%      TUNINGVIEWER, by itself, creates a new TUNINGVIEWER or raises the existing
%      singleton*.
%
%      H = TUNINGVIEWER returns the handle to a new TUNINGVIEWER or the handle to
%      the existing singleton*.
%
%      TUNINGVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TUNINGVIEWER.M with the given input arguments.
%
%      TUNINGVIEWER('Property','Value',...) creates a new TUNINGVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TuningViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TuningViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TuningViewer

% Last Modified by GUIDE v2.5 05-May-2023 11:34:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TuningViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @TuningViewer_OutputFcn, ...
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


% --- Executes just before TuningViewer is made visible.
function TuningViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TuningViewer (see VARARGIN)

% Choose default command line output for TuningViewer
handles.output = hObject;

[Paths] = WhichComputer();
handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
handles.CurrentUnit.Data(1) = NaN;


set(handles.axes10,'Color','none');
set(handles.axes11,'Color','none');
set(handles.axes12,'Color','none');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TuningViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TuningViewer_OutputFcn(hObject, eventdata, handles) 
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
[~, ~, ~, handles.SingleUnits, ~, ...
    ~, handles.SampleRate, ~, ~, ~, ...
    ~, handles.Tuning] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

%% Get all spikes, all units aligned to trial start
[handles.AlignedSpikes, handles.whichtetrode] = TrialAlignedSpikeTimes_Tuning(handles.SingleUnits,handles.Tuning.TTLs);

if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
    % transition markers
    odorTS(1,1) = sum(handles.Tuning.extras.sessionsettings(1,4)); % w.r.t. trial start (motor-settle, pre-odor)
    nLocations = size(handles.Tuning.extras.sequence,2) - 2;
    LocationShifts = 0; 
    for i = 1:nLocations
        if i == 1
            LocationShifts(i,1) = -handles.Tuning.extras.sessionsettings(1,3);
            LocationShifts(i,2) = sum(handles.Tuning.extras.sessionsettings(1,[4,5]));
        else
            LocationShifts(i,1) = LocationShifts(i-1,2);
            LocationShifts(i,2) = LocationShifts(i,1) + sum(handles.Tuning.extras.sessionsettings(1,[4,5]));
        end
    end
    odorTS(1,2) = LocationShifts(end,2);
    handles.TuningTiming.LocationShifts = LocationShifts/1000; % in s
    handles.TuningTiming.Odor = odorTS/1000; % in s
end

handles.NumUnits.String = num2str(size(handles.SingleUnits,2));
handles.CurrentUnit.Data(1) = 1;

UpdatePlots(handles);

% Update handles structure
guidata(hObject, handles);

function UpdatePlots(handles)
whichUnit = handles.CurrentUnit.Data(1);
AlignType = handles.AlignTo.Value;
MyColors1 = brewermap(15,'*PuBu');
MyColors2 = brewermap(15,'*OrRd');
handles.tetrode.String = num2str(handles.whichtetrode(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.whichtetrode(whichUnit,2));
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
    % plot odor ON, location ON periods
    PlotRandomTuningSession(whichUnit, i, handles.AlignedSpikes, handles.TuningTiming, handles.Tuning.extras.sequence, AlignType); %, varargin);
    %[trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, i, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignType, 'plotspikes', 0);
    set(gca, 'XLim', myXlim);
    UpdateUnits(handles);
end

function UpdateUnits(handles)
whichUnit = handles.CurrentUnit.Data(1);
AlignType = handles.AlignTo.Value;
MyColors1 = brewermap(15,'*PuBu');
MyColors2 = brewermap(15,'*OrRd');
handles.tetrode.String = num2str(handles.whichtetrode(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.whichtetrode(whichUnit,2));
switch AlignType
    case {1,2,6}
        myXlim = [-1.2 6];
    case {3,4}
        myXlim = [-5.2 1];
    case 5
        myXlim = [-2.2 5];
end
    
for i = 1:3
    axes(handles.(['axes',num2str(i+9)])); 
    cla reset; 
    set(gca,'color','none');
    hold on
    % plot baseline trials
    [trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, i, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignType, 'plotevents', 0);
    set(gca, 'XLim', myXlim);
    set(gca, 'YLim', handles.(['axes',num2str(i)]).YLim);
    
    if handles.plotPSTH.Value
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
    
    % plot spike amplitudes
    thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
    axes(handles.amplitudeaxes);
    plot(handles.SingleUnits(whichUnit).spikes,thisunitamps,'.');
    hold on
    session_end = handles.TrialInfo.SessionTimestamps(end,2) + handles.TimestampAdjuster;
    line([session_end session_end],get(gca,'YLim'),'Color','k');
    hold off
    % find 
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
UpdateUnits(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in PrevUnit.
function PrevUnit_Callback(hObject, eventdata, handles)
% hObject    handle to PrevUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CurrentUnit.Data(1) = max(handles.CurrentUnit.Data(1)-1,1);
UpdateUnits(handles);
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
UpdateUnits(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in SortReplay.
function SortReplay_Callback(hObject, eventdata, handles)
% hObject    handle to SortReplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SortReplay
UpdatePlots(handles);
UpdateUnits(handles);
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


% --- Executes on button press in plotPSTH.
function plotPSTH_Callback(hObject, eventdata, handles)
% hObject    handle to plotPSTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotPSTH
