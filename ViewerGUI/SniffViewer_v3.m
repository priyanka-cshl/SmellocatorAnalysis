function varargout = SniffViewer_v3(varargin)
% SNIFFVIEWER_V3 MATLAB code for SniffViewer_v3.fig
%      SNIFFVIEWER_V3, by itself, creates a new SNIFFVIEWER_V3 or raises the existing
%      singleton*.
%
%      H = SNIFFVIEWER_V3 returns the handle to a new SNIFFVIEWER_V3 or the handle to
%      the existing singleton*.
%
%      SNIFFVIEWER_V3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNIFFVIEWER_V3.M with the given input arguments.
%
%      SNIFFVIEWER_V3('Property','Value',...) creates a new SNIFFVIEWER_V3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SniffViewer_v3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SniffViewer_v3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SniffViewer_v3

% Last Modified by GUIDE v2.5 04-Mar-2025 14:43:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SniffViewer_v3_OpeningFcn, ...
                   'gui_OutputFcn',  @SniffViewer_v3_OutputFcn, ...
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


% --- Executes just before SniffViewer_v3 is made visible.
function SniffViewer_v3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SniffViewer_v3 (see VARARGIN)

% Choose default command line output for SniffViewer_v3
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

% housekeeping
handles.CurrentUnit.Data(1) = NaN;
handles.SelectedSniffs = [];

for i = 1:4
    currentaxis = handles.(['Spikes',num2str(i)]);
    axes(currentaxis);
    set(currentaxis,'Color','none');
    hold on
    
    handles.(['spikesplot',num2str(i)]) = plot(NaN, NaN, '.k','Markersize', 0.5);
    handles.(['overlayspikesplot',num2str(i)]) = plot(NaN, NaN, '.r','Markersize', 0.5);
end

% Update handles structure
guidata(hObject, handles);

if exist(handles.WhereSession.String)==2
    LoadSession_Callback(hObject, eventdata, handles);
end
% UIWAIT makes SniffViewer_v3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SniffViewer_v3_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
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

tic
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

fprintf('loading data: '); toc

handles = UpdatePlots(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in ReloadUnits.
function ReloadUnits_Callback(hObject, eventdata, handles)
% hObject    handle to ReloadUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% need to reprocess single units 
[SingleUnits] = ReSaveSingleUnits(handles.WhereSession.String);
for whichunit = 1:size(SingleUnits,2)
    handles.AllUnits.Spikes{whichunit}      = SingleUnits(whichunit).spikes; % raw timestamps in OEPS base
    handles.AllUnits.ChannelInfo(whichunit,1:2) = [SingleUnits(whichunit).tetrode SingleUnits(whichunit).id]; % tetrode and phy cluster ID
end
handles.NumUnits.String = num2str(size(handles.AllUnits.Spikes,2));
% Update handles structure
guidata(hObject, handles);

function [handles] = UpdatePlots(hObject, eventdata, handles)

tic

% get sniff time stamps and info for the sniffs we want to plot
[handles.SelectedSniffs] = SelectSniffs(handles.TrialAligned, handles.TrialInfo, [1 2 3], 'includeITI', 1);
                            
% also get air sniffs from tuning
[handles.SelectedSniffs{4}] = SelectSniffs(handles.TuningAligned, handles.TuningInfo, 0, 'includeITI', 1);
                            
% sort the sniffs
for whichodor = 1:3
    handles.SelectedSniffs{whichodor} = ...
        SortSniffs(handles.SelectedSniffs{whichodor}, handles.SortOrder.Value);
end

myXlim = eval(handles.xlims.String);

fprintf('processing data: '); toc

tic
for whichodor = 1:3
    axes(handles.(['Sniffs',num2str(whichodor)]));  %#ok<LAXES>
    cla reset; 
    hold on

    nSniffs = SniffAlignedPlot(handles.SelectedSniffs{whichodor}, [], ...
                'plotspikes', 0, ...
                'alignto', handles.SniffAlignment.Value, ...
                'warptype', handles.WarpType.Value-1);
            
    set(gca, 'XLim', myXlim, 'YLim', [0 nSniffs], 'XTick', [], 'YTick', []);
    
    axes(handles.(['Spikes',num2str(whichodor)]));  %#ok<LAXES>
    set(gca, 'XLim', myXlim, 'YTick', [], 'YLim', handles.(['Sniffs',num2str(whichodor)]).YLim);
            
end
fprintf('plotting sniffs: '); toc

% Update handles structure
guidata(hObject, handles);

UpdateUnits(handles);

function UpdateUnits(handles)

whichUnit = handles.CurrentUnit.Data(1);
handles.tetrode.String = num2str(handles.AllUnits.ChannelInfo(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.AllUnits.ChannelInfo(whichUnit,2));
myXlim = eval(handles.xlims.String);
PSTHOffset = -1000;
FRLims = [0 0];

if handles.MergeUnits.Value
    OverlayUnit = handles.MergeUnit.Data(1);
    handles.tetrode2.String = num2str(handles.AllUnits.ChannelInfo(OverlayUnit,1));
    handles.Cluster2_ID.String = num2str(handles.AllUnits.ChannelInfo(OverlayUnit,2));
end

tic
for whichodor = 1:3
    
    [nSniffs{whichodor},AllFR{whichodor},SpikesPlot{whichodor}] = SniffAlignedPlot(handles.SelectedSniffs{whichodor}, handles.AllUnits.Spikes{whichUnit}, ...
                'plotevents', 0, ...
                'psth', handles.plotPSTH.Value, ...
                'psthoffset', PSTHOffset, ...
                'alignto', handles.SniffAlignment.Value, ...
                'warptype', handles.WarpType.Value-1);

    if handles.MergeUnits.Value
        [~,~,Spikes2Overlay{whichodor}] = SniffAlignedPlot(handles.SelectedSniffs{whichodor}, handles.AllUnits.Spikes{OverlayUnit}, ...
            'plotevents', 0, ...
            'psth', handles.plotPSTH.Value, ...
            'psthoffset', PSTHOffset, ...
            'alignto', handles.SniffAlignment.Value, ...
            'warptype', handles.WarpType.Value-1);
    end

    if handles.plotPSTH.Value
        FRLims(1) = max(FRLims(1), min(cell2mat(cellfun(@min, AllFR{whichodor}, 'UniformOutput', false))));
        FRLims(2) = max(FRLims(2), max(cell2mat(cellfun(@max, AllFR{whichodor}, 'UniformOutput', false))));
    end
end

for whichodor = 1:3
    % plotspikes
    myspikeplot = ['spikesplot',num2str(whichodor)];
    set(handles.(myspikeplot),'XData',SpikesPlot{whichodor}(:,1),'YData',SpikesPlot{whichodor}(:,2));
    overlayspikeplot = ['overlayspikesplot',num2str(whichodor)];
    if handles.MergeUnits.Value
        set(handles.(overlayspikeplot),'XData',Spikes2Overlay{whichodor}(:,1),'YData',Spikes2Overlay{whichodor}(:,2));
    else
        set(handles.(overlayspikeplot),'XData',nan,'YData',nan);
    end
end
drawnow
fprintf('plotting spikes: '); toc

if handles.plotPSTH.Value
    % plot PSTHs
    MyColors1 = brewermap(8,'YlOrRd');
    ColorList(1,:) = Plot_Colors('b'); % ITI sniffs
    ColorList(2,:) = MyColors1(4,:); % snifftype = 0 % approach
    ColorList(3,:) = MyColors1(6,:); % snifftype = 1 % settle
    ColorList(4,:) = MyColors1(7,:); % snifftype = 2 % at target
    
    for whichodor = 1:3
        axes(handles.(['Psth',num2str(whichodor)])); %#ok<LAXES>
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
    %thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
    axes(handles.amplitudeaxes);
    plot(handles.AllUnits(whichUnit).Spikes,handles.AllUnits(whichUnit).SpikeAmps,'.');
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
temp = rem(handles.CurrentUnit.Data(1)+1, str2double(handles.NumUnits.String));
if ~temp
    handles.CurrentUnit.Data(1) = str2double(handles.NumUnits.String);
else
    handles.CurrentUnit.Data(1) = temp;
end
UpdateUnits(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PrevUnit.
function PrevUnit_Callback(hObject, eventdata, handles)
% hObject    handle to PrevUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = rem(handles.CurrentUnit.Data(1)-1, str2double(handles.NumUnits.String));
if ~temp
    handles.CurrentUnit.Data(1) = str2double(handles.NumUnits.String);
else
    handles.CurrentUnit.Data(1) = temp;
end
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

% --- Executes on button press in MergeUnits.
function MergeUnits_Callback(hObject, eventdata, handles)
% hObject    handle to MergeUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
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
        set(handles.(['Sniffs',num2str(i)]).Children, 'Visible', 'on')
    else
        set(handles.(['Sniffs',num2str(i)]).Children, 'Visible', 'off')
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
    set(handles.(['Sniffs',num2str(i)]),'XLim', myXlim);
    set(handles.(['Spikes',num2str(i)]),'XLim', myXlim);
end
% Hints: get(hObject,'String') returns contents of xlims as text
%        str2double(get(hObject,'String')) returns contents of xlims as a double

