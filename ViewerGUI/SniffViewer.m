function varargout = SniffViewer(varargin)
% SNIFFVIEWER MATLAB code for SniffViewer.fig
%      SNIFFVIEWER, by itself, creates a new SNIFFVIEWER or raises the existing
%      singleton*.
%
%      H = SNIFFVIEWER returns the handle to a new SNIFFVIEWER or the handle to
%      the existing singleton*.
%
%      SNIFFVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNIFFVIEWER.M with the given input arguments.
%
%      SNIFFVIEWER('Property','Value',...) creates a new SNIFFVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SniffViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SniffViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SniffViewer

% Last Modified by GUIDE v2.5 31-Aug-2023 14:42:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SniffViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @SniffViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && isempty(strfind(varargin{1},'.mat'))
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SniffViewer is made visible.
function SniffViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SniffViewer (see VARARGIN)

% Choose default command line output for SniffViewer
handles.output = hObject;

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

handles.CurrentUnit.Data(1) = NaN;

set(handles.axes10,'Color','none');
set(handles.axes11,'Color','none');
set(handles.axes12,'Color','none');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SniffViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SniffViewer_OutputFcn(hObject, eventdata, handles) 
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
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop, handles.Tuning] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

handles.SingleUnits = SingleUnits;
handles.TimestampAdjuster = TimestampAdjuster;

%% Get all spikes, all units aligned to trials
[handles.AlignedSniffs, handles.sniffAlignedSpikes, handles.trialAlignedSpikes, handles.whichtetrode] = ...
    SniffAlignedSpikeTimes(SingleUnits,TTLs,size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession);

%% same for replays
if ~isempty(OpenLoop)
    [handles.ReplayAlignedSniffs, handles.SniffAlignedReplaySpikes, handles.ReplayInfo] = ...
        SniffAlignedSpikeTimes_Replays(SingleUnits,TTLs,ReplayTTLs,handles.TrialInfo,OpenLoop,MySession);
else
    handles.ReplayAlignedSniffs = [];
end

%% also for passive tuning
handles.TuningSniffs = PassiveTuningSniffs(handles.Tuning,MySession);

handles.NumUnits.String = num2str(size(SingleUnits,2));
if isnan(handles.CurrentUnit.Data(1)) || handles.CurrentUnit.Data(1)>size(SingleUnits,2)
    handles.CurrentUnit.Data(1) = 1;
end

UpdatePlots(handles);

% Update handles structure
guidata(hObject, handles);

function UpdatePlots(handles)
whichUnit = handles.CurrentUnit.Data(1);
MyColors1 = brewermap(15,'*PuBu');
MyColors2 = brewermap(15,'*OrRd');
handles.tetrode.String = num2str(handles.whichtetrode(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.whichtetrode(whichUnit,2));

myXlim = eval(handles.xlims.String); %[-0.1 1.1];

for i = 1:3
    axes(handles.(['axes',num2str(i)])); 
    cla reset; 
    hold on
    
    % plot baseline trials
    [nSniffs] = PlotSortedSniffs(whichUnit, i, handles.trialAlignedSpikes, handles.AlignedSniffs, ...
                                 handles.TrialInfo, 'plotspikes', 0, 'sortorder', (handles.SortOrder.Value-1), ...
                                 'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1);
                             
    
%     % plot passive replay trials                         
%     if ~isempty(handles.ReplayAlignedSniffs)
%         [nSniffs] = PlotPassiveReplaySniffs(nSniffs, whichUnit, i, handles.SniffAlignedReplaySpikes, handles.ReplayAlignedSniffs, ...
%             handles.ReplayInfo, 'plotspikes', 0, 'sortorder', (handles.SortOrder.Value-1), ...
%             'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1);
%     end
    
%     % add tuning sniffs
%     [nSniffs] = PlotTuningSniffs(whichUnit, i, handles.SingleUnits, handles.TuningSniffs, handles.Tuning.extras.sequence, nSniffs, ...
%         'plotspikes', 0, 'sortorder', (handles.SortOrder.Value-1), ...
%         'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1);
            
    set(gca, 'XLim', myXlim, 'YLim', [0 nSniffs], 'YTick', []);
end
UpdateUnits(handles);

function UpdateUnits(handles)
whichUnit = handles.CurrentUnit.Data(1);
MyColors1 = brewermap(8,'YlOrRd');
MyColors2 = brewermap(15,'*OrRd');
handles.tetrode.String = num2str(handles.whichtetrode(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.whichtetrode(whichUnit,2));

myXlim = eval(handles.xlims.String); %[-0.1 1.1];

myYlims = [];
for i = 1:3
    axes(handles.(['axes',num2str(i+9)])); 
    cla reset; 
    set(gca,'color','none');
    hold on
    % plot baseline trials
    [nSniffs,FR,BinOffset] = PlotSortedSniffs(whichUnit, i, handles.trialAlignedSpikes, handles.AlignedSniffs, ...
                     handles.TrialInfo, 'plotevents', 0, 'sortorder', (handles.SortOrder.Value-1), ...
                     'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1, ...
                     'psth', handles.plotPSTH.Value);
    
%     % plot passive replay trials                         
%     if ~isempty(handles.ReplayAlignedSniffs)
%         [nSniffs] = PlotPassiveReplaySniffs(nSniffs, whichUnit, i, handles.SniffAlignedReplaySpikes, handles.ReplayAlignedSniffs, ...
%             handles.ReplayInfo, 'plotevents', 0, 'sortorder', (handles.SortOrder.Value-1), ...
%             'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1);
%     end
%     % add tuning sniffs
%     [~] = PlotTuningSniffs(whichUnit, i, handles.SingleUnits, handles.TuningSniffs, handles.Tuning.extras.sequence, nSniffs, ...
%         'plotevents', 0, 'sortorder', (handles.SortOrder.Value-1), ...
%         'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1);
    
    set(gca, 'XLim', myXlim, 'YTick', []);
    set(gca, 'YLim', handles.(['axes',num2str(i)]).YLim);
   
    % plot PSTHs
    if handles.plotPSTH.Value
        axes(handles.(['axes',num2str(i+3)]));
        cla reset;
        hold on
        
        for t = 1:size(FR,2)
            %set(groot,'defaultAxesColorOrder',MyColors2);
            if t < 5
                plot((1:size(FR{t},1))*0.002+BinOffset/1000,FR{t},'Linewidth',2,'Color',MyColors1(t+3,:));
            else
                plot((1:size(FR{t},1))*0.002+BinOffset/1000,FR{t},'Linewidth',2,'Color','k');
            end
        end
        set(gca, 'XLim', myXlim);
        myYlims(i,:) = get(gca, 'YLim');
    end
end

if handles.plotPSTH.Value
    for i = 1:3
        set(handles.(['axes',num2str(i+3)]),'YLim', [min(myYlims(:,1)) max(myYlims(:,2))]);
    end
end

% plot spike amplitudes
if handles.spike_amplitudes.Value
    hold off
    thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
    axes(handles.amplitudeaxes);
    plot(handles.SingleUnits(whichUnit).spikes,thisunitamps,'.');
    hold on
    session_end = handles.TrialInfo.SessionTimestamps(end,2) + handles.TimestampAdjuster;
    line([session_end session_end],get(gca,'YLim'),'Color','k');
end

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


% --- Executes on button press in plotPSTH.
function plotPSTH_Callback(hObject, eventdata, handles)
% hObject    handle to plotPSTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotPSTH


% --- Executes on button press in ShadeExhalation.
function ShadeExhalation_Callback(hObject, eventdata, handles)
% hObject    handle to ShadeExhalation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i = 1:3
    if handles.ShadeExhalation.Value
        set(handles.(['axes',num2str(i)]).Children, 'Visible', 'on')
    else
        set(handles.(['axes',num2str(i)]).Children, 'Visible', 'off')
    end
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of ShadeExhalation



function xlims_Callback(hObject, eventdata, handles)
% hObject    handle to xlims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

myXlim = eval(handles.xlims.String); %[-0.1 1.1];

for i = 1:3
    set(handles.(['axes',num2str(i+0)]),'XLim', myXlim);
    set(handles.(['axes',num2str(i+9)]),'XLim', myXlim);
end
% Hints: get(hObject,'String') returns contents of xlims as text
%        str2double(get(hObject,'String')) returns contents of xlims as a double
