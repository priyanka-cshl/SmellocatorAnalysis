function varargout = ProcessSniffTimeStamps_GUI_v3(varargin)
% PROCESSSNIFFTIMESTAMPS_GUI_V3 MATLAB code for ProcessSniffTimeStamps_GUI_v3.fig
%      PROCESSSNIFFTIMESTAMPS_GUI_V3, by itself, creates a new PROCESSSNIFFTIMESTAMPS_GUI_V3 or raises the existing
%      singleton*.
%
%      H = PROCESSSNIFFTIMESTAMPS_GUI_V3 returns the handle to a new PROCESSSNIFFTIMESTAMPS_GUI_V3 or the handle to
%      the existing singleton*.
%
%      PROCESSSNIFFTIMESTAMPS_GUI_V3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCESSSNIFFTIMESTAMPS_GUI_V3.M with the given input arguments.
%
%      PROCESSSNIFFTIMESTAMPS_GUI_V3('Property','Value',...) creates a new PROCESSSNIFFTIMESTAMPS_GUI_V3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProcessSniffTimeStamps_GUI_v3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProcessSniffTimeStamps_GUI_v3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProcessSniffTimeStamps_GUI_v3

% Last Modified by GUIDE v2.5 31-Jul-2025 11:41:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcessSniffTimeStamps_GUI_v3_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcessSniffTimeStamps_GUI_v3_OutputFcn, ...
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


% --- Executes just before ProcessSniffTimeStamps_GUI_v3 is made visible.
function ProcessSniffTimeStamps_GUI_v3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProcessSniffTimeStamps_GUI_v3 (see VARARGIN)

% Choose default command line output for ProcessSniffTimeStamps_GUI_v3
handles.output = hObject;

% some initializations
handles.SniffTrace.Raw          = [];
handles.SniffTrace.Filtered     = [];
handles.SniffTrace.Timestamps   = [];
handles.OdorLocationTrace       = [];
handles.SniffsTS                = [];
handles.SessionLength           = [];
handles.KeepNewSD.BackgroundColor = [0.94 0.94 0.94];
handles.KeepOldSD.BackgroundColor = [0.94 0.94 0.94];
handles.lastdeleted.String = '';
handles.new_detections.Data(:,1) = [];
handles.datamode = 'smellocator';
handles.ProcessedSession = '';

% plots
axes(handles.SniffingRaw);
hold on;
set(gca,'XGrid', 'on');
handles.rawTrace    = plot(nan,nan,'color',Plot_Colors('b'));
handles.peaksRaw    = plot(nan,nan,'ok','MarkerSize',4); 
handles.valleysRaw  = plot(nan,nan,'or','MarkerSize',4);

axes(handles.MFS);
hold on;
set(gca,'YGrid', 'on');
set(gca,'XGrid', 'on');
handles.mfsTrace    = plot(nan,nan,'color',Plot_Colors('pd'));
handles.peaksmfs    = plot(nan,nan,'ok','MarkerSize',4); 
handles.valleysmfs  = plot(nan,nan,'or','MarkerSize',4);

axes(handles.MFS2Therm);
hold on;
set(gca,'XGrid', 'on');
handles.mfs2thermTrace    = plot(nan,nan,'color',Plot_Colors('t'));
handles.peaksmfs2therm    = plot(nan,nan,'ok','MarkerSize',4); 
handles.valleysmfs2therm  = plot(nan,nan,'or','MarkerSize',4);

axes(handles.SniffingFiltered);
hold on;
set(gca,'XGrid', 'on');
handles.filtTrace    = plot(nan,nan,'color',Plot_Colors('b'));
handles.peaksFilt    = plot(nan,nan,'ok','MarkerSize',4); 
handles.valleysFilt  = plot(nan,nan,'or','MarkerSize',4);
handles.peaksNew     = plot(nan,nan,'og','MarkerSize',6); 
handles.valleysNew   = plot(nan,nan,'om','MarkerSize',6); 


% For loading the processed session
[Paths] = WhichComputer();

if ~isempty(varargin)
    if isstruct(varargin{1})
        handles.WhereSession.String = 'direct traces';
        handles.SniffTrace.Timestamps   = varargin{1}.Timestamps{1};
        handles.OdorLocationTrace       = varargin{1}.Motor{1};
        handles.SniffTrace.Raw          = varargin{1}.Sniffs{1};
        handles.ProcessedSession        = varargin{2};
    elseif exist(varargin{1}) == 2
        handles.WhereSession.String = varargin{1};
    else
        MouseName = regexprep(varargin{1},'_(\w+)_processed.mat','');
        handles.WhereSession.String = fullfile(Paths.ProcessedSessions,MouseName,varargin{1});
    end
else
    handles.WhereSession.String = ''; %fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
end

% Update handles structure
guidata(hObject, handles);

LoadSession_Callback(hObject, eventdata, handles);

% UIWAIT makes ProcessSniffTimeStamps_GUI_v3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ProcessSniffTimeStamps_GUI_v3_OutputFcn(hObject, eventdata, handles) 
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

if exist(handles.WhereSession.String)==2
    % get data mode
    if ~isempty(strfind(handles.WhereSession.String,'_processed'))
        handles.datamode = 'smellocator';
    elseif ~isempty(strfind(handles.WhereSession.String,'_cid-processed'))
        handles.datamode = 'cid';
    elseif ~isempty(strfind(handles.WhereSession.String,'quickprocesssniffs'))
        handles.datamode = 'onlyEphys';
        %handles.SDfactor.String = '10';
    else
        disp('unknown data format');
        return;
    end
elseif strcmp(handles.WhereSession.String,'direct traces')
    handles.datamode = 'smellocator_raw';
else
    disp('invalid data file');
    return;
end

switch handles.datamode
    case 'cid'
        % data from the concentration - identity experiments
        % load the sniff traces
        load(handles.WhereSession.String,'SniffData','TTLs');
        handles.SessionLength = SniffData(end,1);
        
        % make equivalent long traces for plotting
        handles.SniffTrace.Timestamps   = SniffData(:,1);
        handles.OdorLocationTrace       = [];
        handles.SniffTrace.Raw          = SniffData(:,3); % filtered MFS trace
        handles.SniffTrace.Raw          = FilterMassFlowSensor(SniffData(:,2));
        handles.SniffTrace.Filtered     = cumsum(handles.SniffTrace.Raw); % expected thermistor trace
        
        % detect sniff timestamps on the expected thermistor trace
        %handles.SDfactor.String = '10';
        [handles.SniffsTS, SniffIndices] = ProcessThermistorData([SniffData(:,1) handles.SniffTrace.Filtered]);
        handles.SniffsTS(:,8:9) = SniffIndices(:,1:2);
        handles.SniffsTS(find(isnan(handles.SniffsTS(:,8))),:) = [];
        handles.SniffsTS(find(handles.SniffsTS(:,8)<=0),:) = [];

    case 'onlyEphys'
        load(handles.WhereSession.String, 'AllSniffs', 'RespirationData'); % AllSniffs: nx13, RespirationData: t x 3(ts, raw, filt, MFS)
        handles.SessionLength = RespirationData(end,1);

        try
            load(fullfile(fileparts(handles.WhereSession.String), 'quickprocessOdorTTLs.mat'));
        catch
            warning('no TTLs found');
        end

        [SniffTimeStamps] = ChunkWiseSniffs(RespirationData(:,[1 3]), 'SDfactor', str2double(handles.SDfactor.String)); % [sniffstart sniffstop nextsniff ~ ~ ~ trialID]
        % remove nans
        SniffTimeStamps(find(isnan(SniffTimeStamps(:,1))),:) = [];

        % remove overlapping sniffs
        handles.SniffsTS = SniffTimeStamps(find(SniffTimeStamps(:,7)>0),:);

        % make equivalent long traces for plotting
        handles.SniffTrace.Timestamps   = RespirationData(:,1);
        handles.OdorLocationTrace       = [];
        handles.SniffTrace.Raw          = RespirationData(:,2); % unfiltered thermistor trace
        handles.SniffTrace.Filtered     = RespirationData(:,3); % filtered thermistor trace
        if size(RespirationData,2) == 4
            handles.SniffTrace.MFS      = 1*(RespirationData(:,4) - 2.5) + 2.5; % filtered MassFlowSensor trace
            handles.SniffTrace.MFS2Therm= cumsum(handles.SniffTrace.MFS-2.5);
            % filter to detrend
            fband = [0.1 30];
            Np    = 4; % filter order
            [b,a] = butter(Np,fband/(500/2)); 
            handles.SniffTrace.MFS2Therm= filtfilt(b,a,handles.SniffTrace.MFS2Therm);

        end
        handles.SniffsTS(find(isnan(handles.SniffsTS(:,8))),:) = [];
        handles.SniffsTS(find(handles.SniffsTS(:,8)<=0),:) = [];

    case 'smellocator'
        % data from the lever task - reprocessing and curating sniffs
        % load the relevant data
        load(handles.WhereSession.String,'Traces','TrialInfo','SampleRate');
        handles.SessionLength = ceil(Traces.Timestamps{end}(end));
        
        % reprocess sniff traces on a trialwise basis
        if isfield(Traces,'Motor')
            Traces.OdorLocation     = Traces.Motor;
        end
        [SniffTimeStamps] = ...
            TrialWiseSniffs(TrialInfo,Traces,'dlgoverride',logical(1)); % [sniffstart sniffstop nextsniff odorlocation sniffslope stimstate trialID]
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

    case 'smellocator_raw'
        % fake a Traces struct
        Traces.Timestamps{1} = handles.SniffTrace.Timestamps;
        Traces.Sniffs{1} = handles.SniffTrace.Raw;
        Traces.OdorLocation{1} = handles.OdorLocationTrace;

        % fake a trialinfo struct
        TrialInfo.Odor = 1;
        TrialInfo.SessionTimestamps = Traces.Timestamps{1}([1 end])';
        TrialInfo.OdorStart = [0 0];
        
        [SniffTimeStamps] = ...
            TrialWiseSniffs(TrialInfo,Traces); % [sniffstart sniffstop nextsniff odorlocation sniffslope stimstate trialID]
        % remove overlapping sniffs
        handles.SniffsTS = SniffTimeStamps(find(SniffTimeStamps(:,end)>0),:);

        % Make long concatenated traces for plotting
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

        handles.SessionLength = TrialInfo.SessionTimestamps(2);
end

handles.SessionDuration.String = num2str(handles.SessionLength);

% update all trace plots

if isfield(handles.SniffTrace,'MFS')
    AxesTag     = {'SniffingFiltered'; 'SniffingRaw'; 'MFS'; 'MFS2Therm'};
else
    AxesTag     = {'SniffingFiltered'; 'SniffingRaw'};
end
PlotTag     = {'filtTrace'; 'rawTrace'; 'mfsTrace'; 'mfs2thermTrace'};
TraceTag    = {'Filtered'; 'Raw'; 'MFS'; 'MFS2Therm'};
PeaksTag    = {'peaksFilt'; 'peaksRaw'; 'peaksmfs'; 'peaksmfs2therm'};
ValleysTag  = {'valleysFilt'; 'valleysRaw'; 'valleysmfs'; 'valleysmfs2therm'};

for i = 1:size(AxesTag,1)
    axes(handles.(AxesTag{i}));
    set(handles.(PlotTag{i}),'XData', handles.SniffTrace.Timestamps, 'YData', handles.SniffTrace.(TraceTag{i}));
    zoom off
    % overlay detected timestamps on the plot
    set(handles.(PeaksTag{i}),'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,8)),...
        'YData',handles.SniffTrace.(TraceTag{i})(handles.SniffsTS(:,8)));
    set(handles.(ValleysTag{i}),'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,9)),...
        'YData',handles.SniffTrace.(TraceTag{i})(handles.SniffsTS(:,9)));

    if i == 1
        set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)],'YTick', []);
        currlims = get(gca,'YLim');

        if exist('TTLs','var')
            % Trial times
            Xvals = [TTLs.Trial(:,1) TTLs.Trial(:,1) nan(size(TTLs.Trial,1),1)]';
            YVals = repmat([2*currlims NaN],size(TTLs.Trial,1),1)';
            plot(Xvals(:),YVals(:),'--','color',Plot_Colors('pd'));

            Xvals = [TTLs.Trial(:,2) TTLs.Trial(:,2) nan(size(TTLs.Trial,1),1)]';
            YVals = repmat([2*currlims NaN],size(TTLs.Trial,1),1)';
            plot(Xvals(:),YVals(:),'-.','color',Plot_Colors('pl'));
        end
        set(gca,'YLim',currlims);
    elseif i == 3
        set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)],'YTick', 2.5,'XTickLabel', {});
    else
        set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)],'YTick', [], 'XTickLabel', {});
    end

end

if strcmp(handles.datamode,'cid')
    line(handles.SniffTrace.Timestamps([1 end]),[0 0],'color','k','Linestyle',':');
end

axes(handles.SniffingFiltered);

% Update handles structure
guidata(hObject, handles);

function [] = UpdatePeakValleyPlots(hObject, eventdata, handles)

% update all plots
TraceTag = {'Filtered'; 'Raw'; 'MFS'; 'MFS2Therm'};
peaksTag = {'Filt'; 'Raw'; 'mfs'; 'mfs2therm'};

% the original points
olddetections = find(handles.SniffsTS(:,8)>=0);

for plotnum = 1:size(TraceTag,1)
    if isfield(handles.SniffTrace, TraceTag{plotnum})
        set(handles.(['peaks',peaksTag{plotnum}]),...
            'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,8)),...
            'YData',handles.SniffTrace.(TraceTag{plotnum})(handles.SniffsTS(olddetections,8)));
        set(handles.(['valleys',peaksTag{plotnum}]),...
            'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,9)),...
            'YData',handles.SniffTrace.(TraceTag{plotnum})(handles.SniffsTS(olddetections,9)));
    end
end

% new detections
if isfield(handles,'SniffsTSnew') & ~isempty(handles.SniffsTSnew)
    newdetections = find(handles.SniffsTSnew(:,9)>0);
    handles.new_detections.Data = round(handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,8)),3,'decimals');
    set(handles.peaksNew, ...
        'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,8)),...
        'YData',handles.SniffTrace.Filtered(handles.SniffsTSnew(newdetections,8)));
    set(handles.valleysNew, ...
        'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,9)),...
        'YData',handles.SniffTrace.Filtered(handles.SniffsTSnew(newdetections,9)));
else
    handles.new_detections.Data = [];
    set(handles.peaksNew, ...
        'XData',[],...
        'YData',[]);
    set(handles.valleysNew, ...
        'XData',[],...
        'YData',[]);
end
guidata(hObject, handles);

% --- Executes on button press in RemovePoints.
function RemovePoints_Callback(hObject, eventdata, handles)
roi = drawrectangle;
SniffTag = 'SniffsTS';

peaks_to_delete = intersect(...
    find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
    find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );

valleys_to_delete = intersect(...
    find(handles.(SniffTag)(:,2) >= roi.Position(1)), ...
    find(handles.(SniffTag)(:,2) <= (roi.Position(1) + roi.Position(3)) ) );

if ~isempty(peaks_to_delete) & numel(peaks_to_delete)==numel(valleys_to_delete)

    if peaks_to_delete(1) <= valleys_to_delete(1)
        handles.SniffsTS(peaks_to_delete(1)-1,3) = handles.SniffsTS(peaks_to_delete(end)+1,1);
        handles.SniffsTS(peaks_to_delete,:) = [];
    else
        handles.SniffsTS(valleys_to_delete(end)+1,1) = handles.SniffsTS(valleys_to_delete(1),1);
        handles.SniffsTS(valleys_to_delete(end)+1,8) = handles.SniffsTS(valleys_to_delete(1),8); % copy indices
        handles.SniffsTS(valleys_to_delete,:) = [];
    end
    guidata(hObject, handles);
    UpdatePeakValleyPlots(hObject, eventdata, handles);
end
delete(roi);
%guidata(hObject, handles);

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
 % Get the key pressed
    pressedKey = eventdata.Key;
    
    % Implement actions based on the key
    switch pressedKey
        case 'd'
            handles.RemovePoints.BackgroundColor = [0.54 0.94 0.54];
            RemovePoints_Callback(hObject, eventdata, handles);
            handles.RemovePoints.BackgroundColor = [0.94 0.94 0.94];
        case 'o'
            % Call a function or perform an action for 'o' key
            disp('Open action triggered!');
            % call open_file_function(handles);
        % Add more cases for other keys
    end

% --- Old version: Executes on button press in RemovePoints.
function RemovePoints_Callback_Old(hObject, eventdata, handles)
% hObject    handle to RemovePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi = drawrectangle;
SniffTag = 'SniffsTSnew';
if ~isfield(handles,SniffTag)
    SniffTag = 'SniffsTS';
end

peaks_to_delete = intersect(...
                    find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
                        find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );

valleys_to_delete = intersect(...
                    find(handles.(SniffTag)(:,2) >= roi.Position(1)), ...
                        find(handles.(SniffTag)(:,2) <= (roi.Position(1) + roi.Position(3)) ) );

if isempty(peaks_to_delete) & isempty(valleys_to_delete) & strcmp(SniffTag,'SniffsTSnew')
    SniffTag = 'SniffsTS';
    peaks_to_delete = intersect(...
        find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
        find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );

    valleys_to_delete = intersect(...
        find(handles.(SniffTag)(:,2) >= roi.Position(1)), ...
        find(handles.(SniffTag)(:,2) <= (roi.Position(1) + roi.Position(3)) ) );
end


if strcmp(SniffTag,'SniffsTS')
    for i = 1:numel(peaks_to_delete)
        handles.SniffsTS(peaks_to_delete(i)-1,3) = handles.SniffsTS(peaks_to_delete(i)+1,1);
        handles.SniffsTS(peaks_to_delete(i),:) = [];
    end
else
    if peaks_to_delete == valleys_to_delete
        handles.lastdeleted.String = mat2str(peaks_to_delete);
        handles.SniffsTSnew(peaks_to_delete,9) = -abs(handles.SniffsTSnew(peaks_to_delete,9));
    elseif peaks_to_delete == valleys_to_delete + 1
        handles.SniffsTSnew(peaks_to_delete,9) = -abs(handles.SniffsTSnew(peaks_to_delete,9));
        handles.SniffsTSnew(valleys_to_delete,9) = -abs(handles.SniffsTSnew(valleys_to_delete,9));
        % unflag the original sniff in the old detections
        [~,m] = min(abs(handles.SniffsTS(:,1)-handles.SniffsTSnew(valleys_to_delete,1)));
        handles.SniffsTS(m,8) = abs(handles.SniffsTS(m,8));

        handles.lastdeleted.String = mat2str(vertcat(valleys_to_delete,peaks_to_delete,-m));
    end
end
guidata(hObject, handles);
UpdatePeakValleyPlots(hObject, eventdata, handles);    
delete(roi);
%guidata(hObject, handles);


% --- Executes on button press in UndoDelete.
function UndoDelete_Callback(hObject, eventdata, handles)
% hObject    handle to UndoDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.lastdeleted.String)
    undelete = eval(handles.lastdeleted.String);
    handles.SniffsTSnew(undelete(undelete>0),9) = abs(handles.SniffsTSnew(undelete(undelete>0),9));
    handles.SniffsTS(abs(undelete(undelete<0)),8) = -abs(handles.SniffsTS(abs(undelete(undelete<0)),8));
    UpdatePeakValleyPlots(hObject, eventdata, handles);
    handles.lastdeleted.String = '';
    guidata(hObject, handles);
end

% --- Executes on button press in ModifyPoint.
function ModifyPoint_Callback(hObject, eventdata, handles)
% hObject    handle to ModifyPoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% user selects which point to edit
roi = drawrectangle;
SniffTag = 'SniffsTSnew';
if ~isfield(handles,SniffTag)
    SniffTag = 'SniffsTS';
end

peaks_to_delete = intersect(...
                    find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
                        find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );
peaks_to_delete = intersect(peaks_to_delete,find(~isnan(handles.(SniffTag)(:,8))));

valleys_to_delete = intersect(...
                    find(handles.(SniffTag)(:,2) >= roi.Position(1)), ...
                        find(handles.(SniffTag)(:,2) <= (roi.Position(1) + roi.Position(3)) ) );
valleys_to_delete = intersect(valleys_to_delete,find(~isnan(handles.(SniffTag)(:,8))));

if isempty([peaks_to_delete; valleys_to_delete])
    % try originally detected sniffs
    SniffTag = 'SniffsTS';
    peaks_to_delete = intersect(...
        find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
        find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );

    valleys_to_delete = intersect(...
        find(handles.(SniffTag)(:,2) >= roi.Position(1)), ...
        find(handles.(SniffTag)(:,2) <= (roi.Position(1) + roi.Position(3)) ) );
end

if ~isempty(peaks_to_delete) & ~isempty(valleys_to_delete)
    warndlg('please select only one peak or valley!','Warning');
    return;
elseif ~isempty(valleys_to_delete)
    ispeak = 0;
    sniff_to_edit = [valleys_to_delete 2];
elseif ~isempty(peaks_to_delete)
    ispeak = 1;
    sniff_to_edit = [peaks_to_delete 1];
end
delete(roi);

% user selects a new point
[x] = ginput(1);
[~,idx] = min(abs(handles.SniffTrace.Timestamps - x(1)));

if ispeak
    if handles.(SniffTag)(sniff_to_edit(1)-1,3) == handles.(SniffTag)(sniff_to_edit(1),1)
        % also edit previous sniff end
        handles.(SniffTag)(sniff_to_edit(1)-1,3) = x(1);
    end
    handles.(SniffTag)(sniff_to_edit(1),1) = x(1);
    handles.(SniffTag)(sniff_to_edit(1),8) = idx;

else
    handles.(SniffTag)(sniff_to_edit(1),2) = x(1);
    handles.(SniffTag)(sniff_to_edit(1),9) = idx;
end

UpdatePeakValleyPlots(hObject, eventdata, handles);
guidata(hObject, handles);

% -----------------------------------------------------------------------

% --- Executes on button press in ZoomON.
function ZoomON_Callback(hObject, eventdata, handles)
axes(handles.SniffingFiltered);
zoom yon
guidata(hObject, handles);

% --- Executes on button press in ZoomOFF.
function ZoomOFF_Callback(hObject, eventdata, handles)
axes(handles.SniffingFiltered);
zoom off
guidata(hObject, handles);

function WindowSize_Callback(hObject, eventdata, handles)
currlims = get(handles.SniffingFiltered,'XLim');
newLims = currlims(1) + [0 str2double(handles.WindowSize.String)];
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);
if isfield(handles.SniffTrace,'MFS')
    set(handles.MFS,'XLim',newLims);
end
if isfield(handles.SniffTrace,'MFS2Therm')
    set(handles.MFS2Therm,'XLim',newLims);
end
guidata(hObject, handles);

% --- Executes on button press in ClearSession.
function ClearSession_Callback(hObject, eventdata, handles)
% hObject    handle to ClearSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on slider movement.
function Scroller_Callback(hObject, eventdata, handles)
newLims = handles.Scroller.Value * handles.SessionLength + ...
    [0 str2double(handles.WindowSize.String)];
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);
if isfield(handles.SniffTrace,'MFS')
    set(handles.MFS,'XLim',newLims);
end
if isfield(handles.SniffTrace,'MFS2Therm')
    set(handles.MFS2Therm,'XLim',newLims);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in NextStretch.
function NextStretch_Callback(hObject, eventdata, handles)
newLims = get(handles.SniffingFiltered,'XLim') + 0.9*str2double(handles.WindowSize.String);
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);
if isfield(handles.SniffTrace,'MFS')
    set(handles.MFS,'XLim',newLims);
end
if isfield(handles.SniffTrace,'MFS2Therm')
    set(handles.MFS2Therm,'XLim',newLims);
end
handles.Scroller.Value = newLims(1)/(handles.SessionLength - str2double(handles.WindowSize.String));
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PreviousStretch.
function PreviousStretch_Callback(hObject, eventdata, handles)
newLims = get(handles.SniffingFiltered,'XLim') - 0.9*str2double(handles.WindowSize.String);
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);
if isfield(handles.SniffTrace,'MFS')
    set(handles.MFS,'XLim',newLims);
end
if isfield(handles.SniffTrace,'MFS2Therm')
    set(handles.MFS2Therm,'XLim',newLims);
end
handles.Scroller.Value = newLims(1)/(handles.SessionLength - str2double(handles.WindowSize.String));
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in RedoStretch.
function RedoStretch_Callback(hObject, eventdata, handles)
roi = drawrectangle;
indices_of_interest = intersect(...
                    find(handles.SniffTrace.Timestamps >= roi.Position(1)), ...
                        find(handles.SniffTrace.Timestamps <= (roi.Position(1) + roi.Position(3)) ) );

redo = 1;
while redo
    redo = redo + 1;
    sdnew = redo*str2double(handles.SDfactor.String);

    stretchTimeStamps = handles.rawTrace.XData(indices_of_interest)';
    stretchThermistor = handles.rawTrace.YData(indices_of_interest)';
    if ~isempty(handles.OdorLocationTrace)
        stretchOdorLocation = handles.OdorLocationTrace(indices_of_interest);
        [thisStretchSniffs, thisStretchIndices] = ...
            ProcessThermistorData([stretchTimeStamps stretchThermistor stretchOdorLocation],'SDfactor',sdnew);
    else
        [thisStretchSniffs, thisStretchIndices] = ...
            ProcessThermistorData([stretchTimeStamps stretchThermistor],'SDfactor',sdnew);
        thisStretchSniffs(find(isnan(thisStretchIndices(:,1))),:) = [];
        if ~isempty(thisStretchSniffs)
            thisStretchIndices(find(isnan(thisStretchIndices(:,1))),:) = [];
        end
    end
    if ~isempty(thisStretchSniffs)
        thisStretchSniffs(:,8:9) = indices_of_interest(thisStretchIndices(:,1:2));
        thisStretchSniffs = thisStretchSniffs(find(thisStretchSniffs(:,end)>0),:);

        % show points on the plot
        set(handles.peaksNew, ...
            'XData',handles.SniffTrace.Timestamps(thisStretchSniffs(:,8)),...
            'YData',handles.SniffTrace.Filtered(thisStretchSniffs(:,8)) );
        set(handles.valleysNew, ...
            'XData',handles.SniffTrace.Timestamps(thisStretchSniffs(:,9)),...
            'YData',handles.SniffTrace.Filtered(thisStretchSniffs(:,9)) );

        % keep the new detections
        answer = questdlg('Keep the new detections?', ...
            'Redo stretch', ...
            'Yes','Redo','Cancel','Yes');

        % Handle response
        switch answer
            case 'Yes'
                % integrate into the current set of sniffs
                whichpeaks = find( (handles.SniffsTS(:,8)>=indices_of_interest(1)) & (handles.SniffsTS(:,8)<=indices_of_interest(end)) );
                whichvalls = find( (handles.SniffsTS(:,9)>=indices_of_interest(1)) & (handles.SniffsTS(:,9)<=indices_of_interest(end)) );
                rows2replace = unique(vertcat(whichpeaks,whichvalls));

                if isempty(rows2replace)
                    % integrate in between sniffs
                    prevsniff = find(handles.SniffsTS(:,1)<thisStretchSniffs(1,1),1,'last');
                    nextsniff = find(handles.SniffsTS(:,1)>thisStretchSniffs(end,1),1,'first');
                    if  nextsniff - prevsniff == 1 && size(thisStretchSniffs,1) == 1
                        chunkA = handles.SniffsTS(1:prevsniff,:);
                        chunkA(end,3) = thisStretchSniffs(1,1);
                        chunkC = handles.SniffsTS(nextsniff:end,:);
                        chunkB = thisStretchSniffs;
                        chunkB(:,4:7) = chunkA(end,4:7);
                        chunkB(1,3) = chunkC(1,1);
                        handles.SniffsTS = [chunkA; chunkB; chunkC];
                    end
                    redo = 0;
                else
                    if numel(rows2replace) > 1
                        % is the last detection before a whole sniff?
                        if thisStretchSniffs(end,2) < handles.SniffsTS(rows2replace(end),1)
                            rows2replace(end,:) = [];
                        end
                    end

                    chunkA = handles.SniffsTS(1:(rows2replace(1)-1),:);
                    chunkC = handles.SniffsTS((rows2replace(end)+1):end,:);
                    chunkB = thisStretchSniffs;
                    chunkB(:,4:7) = repmat(handles.SniffsTS(rows2replace(1),4:7),size(thisStretchSniffs,1),1);
                    chunkB(1,[1 8]) = handles.SniffsTS(rows2replace(1),[1 8]);
                    chunkB(end,3) = chunkC(1,1);
                    handles.SniffsTS = [chunkA; chunkB; chunkC];

                    redo = 0;
                end

            case 'Redo'

            case 'Cancel'
                redo = 0;

        end
        set(handles.peaksNew,'XData',[],'YData',[]);
        set(handles.valleysNew,'XData',[],'YData',[]);
    end
end

UpdatePeakValleyPlots(hObject, eventdata, handles);
delete(roi);
guidata(hObject, handles);

% -----------------------------------------------------------------------

% --- Executes on button press in KeepNewSD.
function KeepNewSD_Callback(hObject, eventdata, handles)
% hObject    handle to KeepNewSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.KeepNewSD.BackgroundColor = [0.94 0.94 0.94];
handles.KeepOldSD.BackgroundColor = [0.94 0.94 0.94];
%uiresume;
handles.SDfactor.String = num2str(2*str2double(handles.SDfactor.String));
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of KeepNewSD

% --- Executes on button press in KeepOldSD.
function KeepOldSD_Callback(hObject, eventdata, handles)
% hObject    handle to KeepOldSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(handles.peaksNew, 'XData',[], 'YData',[]);
% set(handles.valleysNew,  'XData',[], 'YData',[]);
handles.SniffsTSnew = [];
handles.SniffsTS(find(handles.SniffsTS(:,8)<0),8) = -handles.SniffsTS(find(handles.SniffsTS(:,8)<0),8);
UpdatePeakValleyPlots(hObject, eventdata, handles); 
handles.KeepNewSD.BackgroundColor = [0.94 0.94 0.94];
handles.KeepOldSD.BackgroundColor = [0.94 0.94 0.94];
%uiresume;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of KeepOldSD


% --- Executes on button press in SaveSniffs.
function SaveSniffs_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSniffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles] = collatesniffs(handles,1);


% --- Executes when selected cell(s) is changed in new_detections.
function new_detections_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to new_detections (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(eventdata.Indices)
    target_TS = handles.new_detections.Data(eventdata.Indices(1,1),1);%
    newLims = target_TS + [0 str2double(handles.WindowSize.String)] - str2double(handles.WindowSize.String)/2;
    set(handles.SniffingRaw,'XLim',newLims);
    set(handles.SniffingFiltered,'XLim',newLims);
    if isfield(handles.SniffTrace,'MFS')
        set(handles.MFS,'XLim',newLims);
    end
    handles.Scroller.Value = newLims(1)/(handles.SessionLength - str2double(handles.WindowSize.String));
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in FlagStretch.
function FlagStretch_Callback(hObject, eventdata, handles)
% hObject    handle to FlagStretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi = drawrectangle;

SniffTag = 'SniffsTS';

sniffs_to_flag = intersect(...
                    find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
                        find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );
handles.(SniffTag)(sniffs_to_flag,8:9) = -abs(handles.(SniffTag)(sniffs_to_flag,8:9));

UpdatePeakValleyPlots(hObject, eventdata, handles);

% indices_to_delete = intersect(...
%                     find(handles.SniffTrace.Timestamps >= roi.Position(1)), ...
%                         find(handles.SniffTrace.Timestamps <= (roi.Position(1) + roi.Position(3)) ) );
% handles.rawTrace.YData(indices_to_delete) = nan;
                    
% axes(handles.SniffingRaw);
% set(handles.rawTrace,'XData', handles.SniffTrace.Timestamps, 'YData', handles.SniffTrace.Raw);
delete(roi);
guidata(hObject, handles);

% ===================================================================================

% --- Executes on button press in AddSniff.
function AddSniff_Callback(hObject, eventdata, handles)
% hObject    handle to AddSniff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom off
[x,y] = ginput(2); % select a pair of peak and valley

[~,idx1] = min(abs(handles.SniffTrace.Timestamps - x(1)));
[~,idx2] = min(abs(handles.SniffTrace.Timestamps - x(2)));
handles.newSniff(1,[1 8]) = [handles.SniffTrace.Timestamps(idx1) idx1];
handles.newSniff(1,[2 9]) = [handles.SniffTrace.Timestamps(idx2) idx2];

% integrate into the current set of sniffs
whichsniff = find(handles.SniffsTS(:,1)<=handles.newSniff(1,1),1,'last');
handles.newSniff(1,3:7) = handles.SniffsTS(whichsniff,3:7);
handles.SniffsTS(whichsniff,3) = handles.newSniff(1,1);
handles.SniffsTS = [handles.SniffsTS(1:whichsniff,:); handles.newSniff; handles.SniffsTS((whichsniff+1):end,:)];
guidata(hObject, handles);
UpdatePeakValleyPlots(hObject, eventdata, handles); 


% --- Executes on button press in TempSave.
function TempSave_Callback(hObject, eventdata, handles)
SniffDetectionThreshold = str2double(handles.SDfactor.String);
Curated_SniffTimestamps = handles.SniffsTS;
lastTimestamp = floor(handles.SniffingFiltered.XLim(2));
save(handles.WhereSession.String,'SniffDetectionThreshold','Curated_SniffTimestamps','lastTimestamp','-append');
disp(['saved sniffs! @ ',num2str(lastTimestamp), ' seconds']);


% --- Executes on button press in RecoverTemp.
function RecoverTemp_Callback(hObject, eventdata, handles)
% hObject    handle to RecoverTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load(handles.WhereSession.String,'Curated_SniffTimestamps','lastTimestamp');
if ~exist('lastTimestamp','var')
    TSdone = input('Timestamp until which data was curated?');
else
    TSdone = lastTimestamp; 
end
whichidx1 = find(Curated_SniffTimestamps>=TSdone,1,'first');
whichidx2 = find(handles.SniffsTS(:,1)>=TSdone,1,'first'); 
if isequal(Curated_SniffTimestamps(whichidx1,1:3),handles.SniffsTS(whichidx2,1:3))
    handles.SniffsTS = [Curated_SniffTimestamps(1:(whichidx1-1),:); handles.SniffsTS(whichidx2:end,:)]; 
else
    keyboard;
end
guidata(hObject, handles); 
UpdatePeakValleyPlots(hObject, eventdata, handles);
