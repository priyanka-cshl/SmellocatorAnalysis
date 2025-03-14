function varargout = SniffTuningViewer(varargin)
% SNIFFTUNINGVIEWER MATLAB code for SniffTuningViewer.fig
%      SNIFFTUNINGVIEWER, by itself, creates a new SNIFFTUNINGVIEWER or raises the existing
%      singleton*.
%
%      H = SNIFFTUNINGVIEWER returns the handle to a new SNIFFTUNINGVIEWER or the handle to
%      the existing singleton*.
%
%      SNIFFTUNINGVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNIFFTUNINGVIEWER.M with the given input arguments.
%
%      SNIFFTUNINGVIEWER('Property','Value',...) creates a new SNIFFTUNINGVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SniffTuningViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SniffTuningViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SniffTuningViewer

% Last Modified by GUIDE v2.5 14-Mar-2025 14:38:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SniffTuningViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @SniffTuningViewer_OutputFcn, ...
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


% --- Executes just before SniffTuningViewer is made visible.
function SniffTuningViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SniffTuningViewer (see VARARGIN)

% Choose default command line output for SniffTuningViewer
handles.output = hObject;

%% make TABs
%Create tab groupA - motor and odors
handles.tgroupA = uitabgroup('Parent', handles.figure1,'TabLocation', 'top',...
    'Position', [0.125 0.04 0.8500 0.8500] );
handles.tabA1 = uitab('Parent', handles.tgroupA, 'Title', 'Stimulus-wise');
handles.tabA2 = uitab('Parent', handles.tgroupA, 'Title', 'Contanimation');
%Place panels into each tab
set(handles.A1,'Parent',handles.tabA1)
set(handles.A2,'Parent',handles.tabA2)
%Reposition each panel to same location as panel 1
handles.A1.Position = [1.000    0.1   172.000   28.0000];
set(handles.A2,'position',get(handles.A1,'position'));

%% session path
[Paths] = WhichComputer();

if ~isempty(varargin)
    if exist(varargin{1}) == 2
        handles.WhereSession.String = varargin{1};
    else
        handles.WhereSession.String = fullfile(Paths.Wolf.processed,'forWDW',varargin{1});
    end
else
    handles.WhereSession.String = fullfile(Paths.Wolf.processed,'forWDW','Q4_20221112_r0_processed.mat');
end
[~,FileName] = fileparts(handles.WhereSession.String);
handles.MouseName = regexprep(FileName,'_(\w+)_processed.mat','');

%% initializations
for i = 1:5
    currentaxis = handles.(['axes',num2str(i)]);
    axes(currentaxis);
    set(currentaxis,'Color','none','XTick',[],'YTick',[]);
    hold on
    handles.(['spikesplot',num2str(i)]) = plot(NaN, NaN, '.k','Markersize', 0.5);
    handles.(['overlayspikesplot',num2str(i)]) = plot(NaN, NaN, '.r','Markersize', 0.5);
end

% Update handles structure
guidata(hObject, handles);

if exist(handles.WhereSession.String)==2
    LoadSession_Callback(hObject, eventdata, handles);
end

% UIWAIT makes SniffTuningViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SniffTuningViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)

    if isempty(handles.WhereSession.String)
        [Paths] = WhichComputer();
            [WhichSession, SessionPath] = uigetfile(...
                                    fullfile(Paths.Wolf.processed,'forWDW','Q4_20221112_r0_processed.mat'),...
                                    'Select Session');
        handles.WhereSession.String = fullfile(SessionPath,WhichSession);
    end

    % get a list of sniffs, also get units
    [handles.AllSniffs, handles.SniffColumnInfo, SingleUnits] = GetAllSniffs(handles.WhereSession.String);

    for whichunit = 1:size(SingleUnits,2)
        handles.AllUnits.Spikes{whichunit}      = SingleUnits(whichunit).spikes; % raw timestamps in OEPS base
        handles.AllUnits.ChannelInfo(whichunit,1:2) = [SingleUnits(whichunit).tetrode SingleUnits(whichunit).id]; % tetrode and phy cluster ID
    end

    handles.thisUnit.Data(1) = size(handles.AllUnits.Spikes,2);
    if isnan(handles.thisUnit.Data(2)) || handles.thisUnit.Data(2)>handles.thisUnit.Data(1)
        handles.thisUnit.Data(2) = 1;
    end

    handles = UpdatePlots(hObject, eventdata, handles);

    % Update handles structure
    guidata(hObject, handles);


% --- Executes within Load Session and at any unit update call.
function [handles] = UpdatePlots(hObject, eventdata, handles)

    % get sniff time stamps and info for the sniffs we want to plot
    [handles.SelectedSniffs] = ParseSniffsByType(handles.AllSniffs, handles.SortSniffsBy.Value);
    %myXlim = eval(handles.xlims.String);

    handles.maxsniffs = max(cellfun(@length, handles.SelectedSniffs));
    for snifftype = 1:5
        set(handles.(['axes',num2str(snifftype)]), 'YLim',[0 handles.maxsniffs+1]);
    end

    % Update handles structure
    guidata(hObject, handles);

    UpdateUnits(handles);

% --- Executes within Load Session and at any unit update or plot update call.
function UpdateUnits(handles)

    whichUnit = handles.thisUnit.Data(2);
    handles.thisUnit.Data(3) = handles.AllUnits.ChannelInfo(whichUnit,2); % cluster ID
    handles.thisUnit.Data(4) = handles.AllUnits.ChannelInfo(whichUnit,1); % tetrode
    %myXlim = eval(handles.xlims.String);

    if handles.MergeUnits.Value
        OverlayUnit = handles.MergingUnit.Data(1);
        handles.MergingUnit.Data(3) = handles.AllUnits.ChannelInfo(OverlayUnit,1);
        handles.MergingUnit.Data(2) = handles.AllUnits.ChannelInfo(OverlayUnit,2);
    end

    thisUnitSpikes = handles.AllUnits.Spikes{whichUnit};
    
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(handles.SelectedSniffs, thisUnitSpikes);

    for snifftype = 1:5
        myspikeplot = ['spikesplot',num2str(snifftype)];
        SpikesPlot = SpikeRaster{snifftype};
        SpikesPlot(:,2) = abs(SpikesPlot(:,2));
        set(handles.(myspikeplot),'XData',SpikesPlot(:,1),'YData',SpikesPlot(:,2));
        if ~handles.MergeUnits.Value
            set(handles.(['overlayspikesplot',num2str(snifftype)]),'XData',nan,'YData',nan);
        end

    end

    if handles.MergeUnits.Value
        MergeUnitSpikes = handles.AllUnits.Spikes{OverlayUnit};
        [OverlaySpikeRaster, ~] = GetSniffLockedSpikes(handles.SelectedSniffs, MergeUnitSpikes);
        for snifftype = 1:5
            myspikeplot = ['overlayspikesplot',num2str(snifftype)];
            SpikesPlot = OverlaySpikeRaster{snifftype};
            SpikesPlot(:,2) = abs(SpikesPlot(:,2));
            set(handles.(myspikeplot),'XData',SpikesPlot(:,1),'YData',SpikesPlot(:,2));
        end
    end
    drawnow


% --- Executes on button press in NextUnit.
function NextUnit_Callback(hObject, eventdata, handles)
    temp = rem(handles.thisUnit.Data(2)+1, handles.thisUnit.Data(1));
    if ~temp
        handles.thisUnit.Data(2) = handles.thisUnit.Data(1);
    else
        handles.thisUnit.Data(2) = temp;
    end
    UpdateUnits(handles);
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in PrevUnit.
function PrevUnit_Callback(hObject, eventdata, handles)
    temp = rem(handles.thisUnit.Data(2)-1, handles.thisUnit.Data(1));
    if ~temp
        handles.thisUnit.Data(2) = handles.thisUnit.Data(1);
    else
        handles.thisUnit.Data(2) = temp;
    end
    UpdateUnits(handles);
    % Update handles structure
    guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in thisUnit.
function thisUnit_CellEditCallback(hObject, eventdata, handles)
    UpdateUnits(handles);
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on selection change in SortSniffsBy.
function SortSniffsBy_Callback(hObject, eventdata, handles)
    handles = UpdatePlots(hObject, eventdata, handles);
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in MergeUnits.
function MergeUnits_Callback(hObject, eventdata, handles)
    handles = UpdatePlots(hObject, eventdata, handles);
    % Update handles structure
    guidata(hObject, handles);
