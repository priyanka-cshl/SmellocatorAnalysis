function varargout = SniffViewer_v2(varargin)
% SNIFFVIEWER_V2 MATLAB code for SniffViewer_v2.fig
%      SNIFFVIEWER_V2, by itself, creates a new SNIFFVIEWER_V2 or raises the existing
%      singleton*.
%
%      H = SNIFFVIEWER_V2 returns the handle to a new SNIFFVIEWER_V2 or the handle to
%      the existing singleton*.
%
%      SNIFFVIEWER_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNIFFVIEWER_V2.M with the given input arguments.
%
%      SNIFFVIEWER_V2('Property','Value',...) creates a new SNIFFVIEWER_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SniffViewer_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SniffViewer_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SniffViewer_v2

% Last Modified by GUIDE v2.5 15-Feb-2024 06:04:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SniffViewer_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @SniffViewer_v2_OutputFcn, ...
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


% --- Executes just before SniffViewer_v2 is made visible.
function SniffViewer_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SniffViewer_v2 (see VARARGIN)

% Choose default command line output for SniffViewer_v2
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

handles.SelectedSniffs = [];

set(handles.axes10,'Color','none');
set(handles.axes11,'Color','none');
set(handles.axes12,'Color','none');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SniffViewer_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SniffViewer_v2_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles) %#ok<DEFNU>
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
[handles.TrialAligned, handles.TrialInfo, ...
    handles.ReplayAligned, handles.ReplayInfo, ...
    handles.TuningAligned, handles.TuningInfo, ...
    handles.AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

handles.NumUnits.String = num2str(size(handles.AllUnits.Spikes,2));
if isnan(handles.CurrentUnit.Data(1)) || handles.CurrentUnit.Data(1)>size(handles.AllUnits.Spikes,2)
    handles.CurrentUnit.Data(1) = 1;
end

handles = UpdatePlots(handles);

% Update handles structure
guidata(hObject, handles);

function [handles] = UpdatePlots(handles)

% get sniff time stamps and info for the sniffs we want to plot
[handles.SelectedSniffs] = SelectSniffs(handles.TrialAligned, handles.TrialInfo, [1 2 3], ...
                                'includeITI', 1);
                            
% sort the sniffs
for whichodor = 1:3
    handles.SelectedSniffs{whichodor} = ...
        SortSniffs(handles.SelectedSniffs{whichodor}, handles.SortOrder.Value);
end

myXlim = eval(handles.xlims.String);

for whichodor = 1:3
    axes(handles.(['axes',num2str(whichodor)]));  %#ok<LAXES>
    cla reset; 
    hold on
    
    nSniffs = SniffAlignedPlot(handles.SelectedSniffs{whichodor}, [], ...
                'plotspikes', 0, ...
                'alignto', handles.SniffAlignment.Value, ...
                'warptype', handles.WarpType.Value-1);
            
    set(gca, 'XLim', myXlim, 'YLim', [0 nSniffs], 'XTick', [], 'YTick', []);
            
end

UpdateUnits(handles);

function UpdateUnits(handles)
whichUnit = handles.CurrentUnit.Data(1);
handles.tetrode.String = num2str(handles.AllUnits.ChannelInfo(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.AllUnits.ChannelInfo(whichUnit,2));
myXlim = eval(handles.xlims.String);
PSTHOffset = -1000;
FRLims = [0 0];

for whichodor = 1:3
    axes(handles.(['axes',num2str(whichodor+9)]));  %#ok<LAXES>
    cla reset; 
    set(gca,'color','none');
    hold on
    
    [nSniffs,AllFR{whichodor}] = SniffAlignedPlot(handles.SelectedSniffs{whichodor}, handles.AllUnits.Spikes{whichUnit}, ...
                'plotevents', 0, ...
                'psth', handles.plotPSTH.Value, ...
                'psthoffset', PSTHOffset, ...
                'alignto', handles.SniffAlignment.Value, ...
                'warptype', handles.WarpType.Value-1);
    
    set(gca, 'XLim', myXlim, 'YTick', []);
    set(gca, 'YLim', handles.(['axes',num2str(whichodor)]).YLim);
    
    if handles.plotPSTH.Value
        FRLims(1) = max(FRLims(1), min(cell2mat(cellfun(@min, AllFR{whichodor}, 'UniformOutput', false))));
        FRLims(2) = max(FRLims(2), max(cell2mat(cellfun(@max, AllFR{whichodor}, 'UniformOutput', false))));
    end
end

% plot PSTHs
MyColors1 = brewermap(8,'YlOrRd');
ColorList(1,:) = Plot_Colors('b'); % ITI sniffs
ColorList(2,:) = MyColors1(4,:); % snifftype = 0 % approach
ColorList(3,:) = MyColors1(6,:); % snifftype = 1 % settle
ColorList(4,:) = MyColors1(7,:); % snifftype = 2 % at target

if handles.plotPSTH.Value
    for whichodor = 1:3
        axes(handles.(['axes',num2str(whichodor+3)])); %#ok<LAXES>
        cla reset;
        hold on
        FR = AllFR{whichodor};
        
        for t = 1:size(FR,2)
            if t < 5
                plot((1:size(FR{t},1))*0.002+PSTHOffset/1000,FR{t},'Linewidth',2,'Color',ColorList(t,:));
            else
                plot((1:size(FR{t},1))*0.002+PSTHOffset/1000,FR{t},'Linewidth',2,'Color','k');
            end
        end
        set(gca, 'XLim', myXlim, 'YLim', FRLims);
        if whichodor>1
            set(gca,'YTickLabel', {});
        end
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
