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
[~, ~, handles.TrialInfo, handles.SingleUnits, ~, ...
    ~, handles.SampleRate, handles.TimestampAdjuster, ~, ~, ...
    ~, handles.Tuning] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

%% Get all spikes, all units aligned to trial start
[handles.AlignedSpikes, handles.whichtetrode] = TrialAlignedSpikeTimes_Tuning(handles.SingleUnits,handles.Tuning.TTLs);

if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
    % transition markers
    odorTS(1,1) = handles.Tuning.extras.sessionsettings(1,4); % w.r.t. trial start (pre-odor)
    nLocations = size(handles.Tuning.extras.sequence,2) - 2;
    LocationShifts = 0; 
    for i = 1:nLocations
        if i == 1
            LocationShifts(i,1) = -handles.Tuning.extras.sessionsettings(1,3); % w.r.t. trial start (settle)
            LocationShifts(i,2) = sum(handles.Tuning.extras.sessionsettings(1,[4,5])); % w.r.t. trial start (pre + odor)
        else
            LocationShifts(i,1) = LocationShifts(i-1,2);
            LocationShifts(i,2) = LocationShifts(i,1) + sum(handles.Tuning.extras.sessionsettings(1,[3,5])); % settle + odor
        end
    end
    odorTS(1,2) = LocationShifts(end,2);
    handles.TuningTiming.LocationShifts = LocationShifts/1000; % in s
    handles.TuningTiming.Odor = odorTS/1000; % in s
end

handles.NumUnits.String = num2str(size(handles.SingleUnits,2));
handles.CurrentUnit.Data(1) = 1;

% Update AlignTo List
AllLocations = unique(handles.Tuning.extras.sequence(:,3:end));
AllLocations(isnan(AllLocations)) = [];
for i = 1:numel(AllLocations)
    handles.AlignTo.String{i+2,1} = num2str(AllLocations(i)); 
end

UpdatePlots(handles);

% Update handles structure
guidata(hObject, handles);

function UpdatePlots(handles)
whichUnit = handles.CurrentUnit.Data(1);
AlignType = handles.AlignTo.Value;
if AlignType>2
    AlignType = 1000+str2double(handles.AlignTo.String{AlignType,1});
end
LocationDuration = mode(diff(handles.TuningTiming.LocationShifts'));
HalfLocations = size(handles.TuningTiming.LocationShifts,1)/2;
switch AlignType
    case 1 % Trial On
        % need an extra offset 
        if handles.TuningTiming.LocationShifts(1,2) > diff(handles.TuningTiming.LocationShifts(2,:))
            my_start_offset = handles.TuningTiming.LocationShifts(1,2) - diff(handles.TuningTiming.LocationShifts(2,:));
            % normalize my location duration to correctly shift the image
            % of location plots
            my_start_offset = my_start_offset/LocationDuration;
            myXlim(1) = 0.5 + my_start_offset;
            myXlim(2) = 0.5 + my_start_offset + HalfLocations*2*LocationDuration;
        end
        
%         myXlim = handles.TuningTiming.LocationShifts([1,end],2)';
%         myXlim(1) = myXlim(1) - LocationDuration;
    case 2 % Odor On
        % first calculate extra offset w.r.t. Trial Start 
        if handles.TuningTiming.LocationShifts(1,2) > diff(handles.TuningTiming.LocationShifts(2,:))
            my_start_offset = handles.TuningTiming.LocationShifts(1,2) - diff(handles.TuningTiming.LocationShifts(2,:));
            % normalize my location duration to correctly shift the image
            % of location plots
            my_start_offset = my_start_offset/LocationDuration;
            myXlim(1) = 0.5 + my_start_offset;
            myXlim(2) = 0.5 + my_start_offset + HalfLocations*2*LocationDuration;
        end
    otherwise
        myXlim =  [0 HalfLocations*2] + 0.5;
        %myXlim =  [0 HalfLocations*2*LocationDuration];
end
for i = 1:4
    axes(handles.(['axes',num2str(i)])); 
    cla reset; 
    PlotRandomTuningSession(whichUnit, i, handles.AlignedSpikes, handles.TuningTiming, handles.Tuning.extras.sequence, AlignType, myXlim, 'plotspikes', 0); %, varargin);
    colormap(handles.(['axes',num2str(i)]),brewermap([500],'rdbu'));
    set(gca, 'XLim', myXlim, 'XTick', [], 'YTick', []);
    UpdateUnits(handles);
end

function UpdateUnits(handles)
whichUnit = handles.CurrentUnit.Data(1);
AlignType = handles.AlignTo.Value;
if AlignType>2
    AlignType = 1000+str2double(handles.AlignTo.String{AlignType,1});
end
handles.tetrode.String = num2str(handles.whichtetrode(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.whichtetrode(whichUnit,2));
LocationDuration = mode(diff(handles.TuningTiming.LocationShifts'));
HalfLocations = size(handles.TuningTiming.LocationShifts,1)/2;
switch AlignType
    case {1,2}
        myXlim = handles.TuningTiming.LocationShifts([1,end],2)';
        myXlim(1) = myXlim(1) - LocationDuration;
    otherwise
        myXlim = LocationDuration*[-floor(HalfLocations) ceil(HalfLocations)];
end
   
PSTHLims = [];
for i = 1:4
    axes(handles.(['axes',num2str(i+9)])); 
    cla reset; 
    set(gca,'color','none');
    hold on
    [FRs, BinOffset] = PlotRandomTuningSession(whichUnit, i, handles.AlignedSpikes, handles.TuningTiming, handles.Tuning.extras.sequence, AlignType, myXlim, 'plotevents', 0); %, varargin);
    myYlim = get(handles.(['axes',num2str(i)]),'YLim') - 0.5;
    set(gca, 'XLim', myXlim, 'XTick', [], 'YLim', myYlim, 'YTick', [], 'TickDir', 'out');
    
    if handles.plotPSTH.Value
        axes(handles.(['axes',num2str(i+4)]));
        cla reset;
        hold on
        plot((1:size(FRs,1))*0.002+BinOffset/1000,FRs,'Color','k','Linewidth',1);
        
        set(gca, 'XLim', myXlim);
        PSTHLims(i, :) = get(gca, 'YLim');
    end
end

if handles.plotPSTH.Value
    myYlim = [min(PSTHLims(:,1)) max(PSTHLims(:,2))];
    for i = 1:4
        axes(handles.(['axes',num2str(i+4)]));
        set(gca, 'YLim', myYlim, 'TickDir', 'out', 'TickLength',[0.002 0.025]);
    end
end

% plot spike amplitudes
thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
axes(handles.amplitudeaxes);
hold off
plot(handles.SingleUnits(whichUnit).spikes,thisunitamps,'.b');
hold on
session_end = handles.TrialInfo.SessionTimestamps(end,2) + handles.TimestampAdjuster;
line([session_end session_end],get(gca,'YLim'),'Color','k');


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
