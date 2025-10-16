function varargout = ProcessSniffTimeStamps_GUI_MFS(varargin)
% PROCESSSNIFFTIMESTAMPS_GUI_MFS MATLAB code for ProcessSniffTimeStamps_GUI_MFS.fig
%      PROCESSSNIFFTIMESTAMPS_GUI_MFS, by itself, creates a new PROCESSSNIFFTIMESTAMPS_GUI_MFS or raises the existing
%      singleton*.
%
%      H = PROCESSSNIFFTIMESTAMPS_GUI_MFS returns the handle to a new PROCESSSNIFFTIMESTAMPS_GUI_MFS or the handle to
%      the existing singleton*.
%
%      PROCESSSNIFFTIMESTAMPS_GUI_MFS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCESSSNIFFTIMESTAMPS_GUI_MFS.M with the given input arguments.
%
%      PROCESSSNIFFTIMESTAMPS_GUI_MFS('Property','Value',...) creates a new PROCESSSNIFFTIMESTAMPS_GUI_MFS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProcessSniffTimeStamps_GUI_MFS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProcessSniffTimeStamps_GUI_MFS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProcessSniffTimeStamps_GUI_MFS

% Last Modified by GUIDE v2.5 06-Aug-2025 13:11:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcessSniffTimeStamps_GUI_MFS_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcessSniffTimeStamps_GUI_MFS_OutputFcn, ...
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


% --- Executes just before ProcessSniffTimeStamps_GUI_MFS is made visible.
function ProcessSniffTimeStamps_GUI_MFS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProcessSniffTimeStamps_GUI_MFS (see VARARGIN)

% Choose default command line output for ProcessSniffTimeStamps_GUI_MFS
handles.output = hObject;

% some initializations
if size(varargin,2) == 2
    switch varargin{2}
        case 1
            handles.primarySniffTS          = 'SniffsMFS'; % 'SniffsTS'
            handles.primaryAxes             = 'MFS2Therm';
            handles.primaryRawSniffTrace    = 'mfs2thermTrace';
            handles.primarySniffTrace       = 'MFS2Therm';
        case 2
            handles.primarySniffTS          = 'SniffsTS'; % 'SniffsTS'
            handles.primaryAxes             = 'SniffingFiltered';
            handles.primaryRawSniffTrace    = 'rawTrace';
            handles.primarySniffTrace       = 'Filtered';
    end
else
    handles.primarySniffTS          = 'SniffsMFS'; % 'SniffsTS'
    handles.primaryAxes             = 'MFS2Therm';
    handles.primaryRawSniffTrace    = 'mfs2thermTrace';
    handles.primarySniffTrace       = 'MFS2Therm';
end

handles.SniffTrace.Raw          = [];
handles.SniffTrace.Filtered     = [];
handles.SniffTrace.Timestamps   = [];
handles.OdorLocationTrace       = [];
handles.SniffsTS                = [];
handles.SniffsMFS               = [];
handles.SessionLength           = [];
handles.KeepNewSD.BackgroundColor = [0.94 0.94 0.94];
handles.KeepOldSD.BackgroundColor = [0.94 0.94 0.94];
handles.lastdeleted.String = '';
handles.datamode = 'smellocator';
handles.ProcessedSession = '';
handles.axesPositions = [];

% plots
axes(handles.SniffingRaw);
hold on;
set(gca,'XGrid', 'on');
handles.rawTrace    = plot(nan,nan,'color',Plot_Colors('b'));
handles.peaksRaw    = plot(nan,nan,'ok','MarkerSize',4); 
handles.valleysRaw  = plot(nan,nan,'or','MarkerSize',4);
handles.axesPositions(1,:) = get(gca,'Position');

axes(handles.MFS);
hold on;
set(gca,'YGrid', 'on');
set(gca,'XGrid', 'on');
handles.mfsTrace    = plot(nan,nan,'color',Plot_Colors('t'));
handles.peaksmfs    = plot(nan,nan,'o','MarkerSize',4,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','none'); 
handles.valleysmfs  = plot(nan,nan,'o','MarkerSize',4,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','none');
handles.zeromfs     = plot(nan,nan,'x','Color',Plot_Colors('p'),'MarkerSize',8,'LineWidth',1); 
handles.pairedTherm = plot(nan,nan,'x','Color',Plot_Colors('b'),'MarkerSize',8,'LineWidth',1); 
handles.peaksBmfs   = plot(nan,nan,'ok','MarkerSize',4,'LineWidth',1.5); 
handles.valleysBmfs = plot(nan,nan,'or','MarkerSize',4,'LineWidth',1.5);
handles.pointercarat= plot(nan,nan,'vk'); 
handles.axesPositions(3,:) = get(gca,'Position');

axes(handles.MFS2Therm);
hold on;
set(gca,'XGrid', 'on');
handles.mfs2thermTrace    = plot(nan,nan,'color',Plot_Colors('t'));
handles.peaksmfs2therm    = plot(nan,nan,'o','MarkerSize',4,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','none'); 
handles.valleysmfs2therm  = plot(nan,nan,'o','MarkerSize',4,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','none');
handles.peaksBmfs2therm    = plot(nan,nan,'ok','MarkerSize',4,'LineWidth',1.5); 
handles.valleysBmfs2therm  = plot(nan,nan,'or','MarkerSize',4,'LineWidth',1.5);
if strcmp(handles.primarySniffTS,'SniffsMFS')
    handles.peaksNew     = plot(nan,nan,'og','MarkerSize',6); 
    handles.valleysNew   = plot(nan,nan,'om','MarkerSize',6); 
end
handles.axesPositions(4,:) = get(gca,'Position');

axes(handles.SniffingFiltered);
hold on;
set(gca,'XGrid', 'on');
handles.filtTrace    = plot(nan,nan,'color',Plot_Colors('b'));
handles.peaksFilt    = plot(nan,nan,'ok','MarkerSize',4); 
handles.valleysFilt  = plot(nan,nan,'or','MarkerSize',4);
handles.pointerline  = line([nan nan],[nan nan],'Color','k'); 
if strcmp(handles.primarySniffTS,'SniffsTS')
    handles.peaksNew     = plot(nan,nan,'og','MarkerSize',6); 
    handles.valleysNew   = plot(nan,nan,'om','MarkerSize',6); 
end
handles.axesPositions(2,:) = get(gca,'Position');

% reorder axes if needed
if strcmp(handles.primarySniffTS,'SniffsTS')
    AxesTag     = {'MFS'; 'MFS2Therm'; 'SniffingRaw'; 'SniffingFiltered'};
    for i = 1:numel(AxesTag)
        set(handles.(AxesTag{i}),'Position',handles.axesPositions(i,:));
    end
end

% For loading the processed session
[Paths] = WhichComputer();
if ~isempty(varargin)
        handles.WhereSession.String = varargin{1};
else
    handles.WhereSession.String = ''; %fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
end

% Update handles structure
guidata(hObject, handles);
LoadSession_Callback(hObject, eventdata, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ProcessSniffTimeStamps_GUI_MFS_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)
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
        handles.SDfactor.String = '2.5';
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
    case 'onlyEphys'
        load(handles.WhereSession.String, 'AllSniffs', 'RespirationData'); % AllSniffs: nx13, RespirationData: t x 3(ts, raw, filt, MFS)
        handles.SessionLength = RespirationData(end,1);

        try
            load(fullfile(fileparts(handles.WhereSession.String), 'quickprocessOdorTTLs.mat'));
        catch
            warning('no TTLs found');
        end

        try 
            load(handles.WhereSession.String,'CuratedSniffTimestamps');
            handles.SniffsTS = CuratedSniffTimestamps;
        catch
            [SniffTimeStamps] = ChunkWiseSniffs(RespirationData(:,[1 3]), 'SDfactor', str2double(handles.SDfactor.String)); % [sniffstart sniffstop nextsniff ~ ~ ~ trialID]
            % remove nans
            SniffTimeStamps(find(isnan(SniffTimeStamps(:,1))),:) = [];

            % remove overlapping sniffs
            handles.SniffsTS = SniffTimeStamps(find(SniffTimeStamps(:,7)>0),:);
        end
        handles.SniffsTS(find(isnan(handles.SniffsTS(:,8))),:) = [];
        handles.SniffsTS(find(handles.SniffsTS(:,8)<=0),:) = [];

        % make equivalent long traces for plotting
        handles.SniffTrace.Timestamps   = RespirationData(:,1);
        handles.OdorLocationTrace       = [];
        handles.SniffTrace.Raw          = RespirationData(:,2); % unfiltered thermistor trace
        handles.SniffTrace.Filtered     = RespirationData(:,3); % filtered thermistor trace

        if size(RespirationData,2) == 4
            handles.SniffTrace.MFS      = RespirationData(:,4) - 2.5; % unfiltered MassFlowSensor trace
            % filter
            fband = [0.1 30];
            Np    = 4; % filter order
            [b,a] = butter(Np,fband/(500/2)); 
            
            % predicted thermistor trace
            %handles.SniffTrace.MFS2Therm= cumsum(filtfilt(b,a,handles.SniffTrace.MFS));
            handles.SniffTrace.MFS2Therm= filtfilt(b,a,cumsum(handles.SniffTrace.MFS));

            % detect peaks and valleys on the predicted thermistor trace
            [SniffTimeStamps] = ChunkWiseSniffs([RespirationData(:,1) handles.SniffTrace.MFS2Therm], 'SDfactor', str2double(handles.SDfactor.String)); % [sniffstart sniffstop nextsniff ~ ~ ~ trialID]
            % remove nans
            SniffTimeStamps(find(isnan(SniffTimeStamps(:,1))),:) = [];

            % remove overlapping sniffs
            handles.SniffsMFS = SniffTimeStamps(find(SniffTimeStamps(:,7)>0),:);

            handles.SniffsMFS(find(isnan(handles.SniffsMFS(:,8))),:) = [];
            handles.SniffsMFS(find(handles.SniffsMFS(:,8)<=0),:) = [];

            handles.SniffsMFS(:,4:7) = nan;
        end
        
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
PeaksTag    = {'peaksFilt'; 'peaksRaw'; 'peaksmfs'; 'peaksmfs2therm'; 'peaksBmfs'; 'peaksBmfs2therm'};
ValleysTag  = {'valleysFilt'; 'valleysRaw'; 'valleysmfs'; 'valleysmfs2therm' ;'valleysBmfs'; 'valleysBmfs2therm'};

for i = 1:size(AxesTag,1)
    axes(handles.(AxesTag{i}));
    set(handles.(PlotTag{i}),'XData', handles.SniffTrace.Timestamps, 'YData', handles.SniffTrace.(TraceTag{i}));
    zoom off
    % overlay detected timestamps on the plot - from the thermistor
    set(handles.(PeaksTag{i}),'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,8)),...
        'YData',handles.SniffTrace.(TraceTag{i})(handles.SniffsTS(:,8)));
    set(handles.(ValleysTag{i}),'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,9)),...
        'YData',handles.SniffTrace.(TraceTag{i})(handles.SniffsTS(:,9)));

    if i > 2
        % overlay MFS detections
        set(handles.(PeaksTag{i+2}),'XData',handles.SniffTrace.Timestamps(handles.SniffsMFS(:,8)),...
            'YData',handles.SniffTrace.(TraceTag{i})(handles.SniffsMFS(:,8)));
        set(handles.(ValleysTag{i+2}),'XData',handles.SniffTrace.Timestamps(handles.SniffsMFS(:,9)),...
            'YData',handles.SniffTrace.(TraceTag{i})(handles.SniffsMFS(:,9)));
    end
    if strcmp(handles.primarySniffTS,'SniffsMFS')
        lowest_axes = 4;
    else
        lowest_axes = 1;
    end
    % axes limits and labels
    if i == lowest_axes
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
        set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)],'YTick', 0,'XTickLabel', {});
        set(handles.(PeaksTag{i}),'XData',[],'YData',[]);
        set(handles.(ValleysTag{i}),'XData',[],'YData',[]);
    else
        set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)],'YTick', [], 'XTickLabel', {});
    end

end

% ste primary working axes
axes(handles.(handles.primaryAxes));

% Update handles structure
guidata(hObject, handles);

function [] = UpdatePeakValleyPlots(hObject, eventdata, handles)

% update all plots
TraceTag = {'Filtered'; 'Raw'; 'MFS'; 'MFS2Therm'};
peaksTag = {'Filt'; 'Raw'; 'mfs'; 'mfs2therm'; 'Bmfs'; 'Bmfs2therm'};
SniffTag = handles.primarySniffTS;

% the original points
olddetections = find(handles.(SniffTag)(:,8)>=0);

if strcmp(handles.primarySniffTS,'SniffsMFS')
    for plotnum = 3:size(TraceTag,1)
        if isfield(handles.SniffTrace, TraceTag{plotnum}) & plotnum >= 3
            set(handles.(['peaks',peaksTag{plotnum+2}]),...
                'XData',handles.SniffTrace.Timestamps(handles.(SniffTag)(olddetections,8)),...
                'YData',handles.SniffTrace.(TraceTag{plotnum})(handles.(SniffTag)(olddetections,8)));
            set(handles.(['valleys',peaksTag{plotnum+2}]),...
                'XData',handles.SniffTrace.Timestamps(handles.(SniffTag)(olddetections,9)),...
                'YData',handles.SniffTrace.(TraceTag{plotnum})(handles.(SniffTag)(olddetections,9)));
        end
        if plotnum == 3
            set(handles.zeromfs,'XData',handles.SniffsMFS(:,4),'YData',0*handles.SniffsMFS(:,4));
            set(handles.pairedTherm,'XData',handles.SniffsMFS(:,5),'YData',handles.SniffsMFS(:,6));
        end
    end
else
    for plotnum = 1:4
        if isfield(handles.SniffTrace, TraceTag{plotnum})
            set(handles.(['peaks',peaksTag{plotnum}]),...
                'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,8)),...
                'YData',handles.SniffTrace.(TraceTag{plotnum})(handles.SniffsTS(olddetections,8)));
            set(handles.(['valleys',peaksTag{plotnum}]),...
                'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,9)),...
                'YData',handles.SniffTrace.(TraceTag{plotnum})(handles.SniffsTS(olddetections,9)));
        end
    end
end

% new detections
if isfield(handles,'SniffsTSnew') & ~isempty(handles.SniffsTSnew)
    newdetections = find(handles.SniffsTSnew(:,9)>0);
    handles.new_detections.Data = round(handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,8)),3,'decimals');
    set(handles.peaksNew, ...
        'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,8)),...
        'YData',handles.SniffTrace.(handles.primarySniffTrace)(handles.SniffsTSnew(newdetections,8)));
    set(handles.valleysNew, ...
        'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,9)),...
        'YData',handles.SniffTrace.(handles.primarySniffTrace)(handles.SniffsTSnew(newdetections,9)));
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
SniffTag = handles.primarySniffTS;

peaks_to_delete = intersect(...
    find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
    find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );

valleys_to_delete = intersect(...
    find(handles.(SniffTag)(:,2) >= roi.Position(1)), ...
    find(handles.(SniffTag)(:,2) <= (roi.Position(1) + roi.Position(3)) ) );

if ~isempty(peaks_to_delete) & numel(peaks_to_delete)==numel(valleys_to_delete)

    if peaks_to_delete(1) <= valleys_to_delete(1)
        handles.(SniffTag)(peaks_to_delete(1)-1,3) = handles.(SniffTag)(peaks_to_delete(end)+1,1);
        handles.(SniffTag)(peaks_to_delete,:) = [];
    else
        handles.(SniffTag)(valleys_to_delete(end)+1,1) = handles.(SniffTag)(valleys_to_delete(1),1);
        handles.(SniffTag)(valleys_to_delete(end)+1,8) = handles.(SniffTag)(valleys_to_delete(1),8); % copy indices
        handles.(SniffTag)(valleys_to_delete,:) = [];
    end
    guidata(hObject, handles);
    UpdatePeakValleyPlots(hObject, eventdata, handles);
end
delete(roi);
%guidata(hObject, handles);

% --- Executes on button press in ModifyPoint.
function ModifyPoint_Callback(hObject, eventdata, handles)
% user selects which point to edit
roi = drawrectangle;
SniffTag = 'SniffsTSnew';
if ~isfield(handles,SniffTag)
    SniffTag = handles.primarySniffTS;
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
    SniffTag = handles.primarySniffTS;
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
% --- Executes on button press in ModifyZeroCross.
function ModifyZeroCross_Callback(hObject, eventdata, handles)
% hObject    handle to ModifyZeroCross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

function GetZeroCrossings_ButtonDownFcn(hObject, eventdata, handles)
handles.roiPosition = get(handles.SniffingFiltered,'XLim');
GetZeroCrossings_Callback(hObject, eventdata, handles)

% --- Executes on button press in GetZeroCrossings.
function GetZeroCrossings_Callback(hObject, eventdata, handles)
axes(handles.MFS);
if ~isfield(handles,'roiPosition')
    roi = drawrectangle;
    sniffs_of_interest = intersect(...
                        find(handles.SniffsMFS(:,1) >= roi.Position(1)), ...
                            find(handles.SniffsMFS(:,2) <= (roi.Position(1) + roi.Position(3)) ) );
else
    sniffs_of_interest = intersect(...
                        find(handles.SniffsMFS(:,1) >= handles.roiPosition(1)), ...
                            find(handles.SniffsMFS(:,2) <= (handles.roiPosition(2)) ) );
end

mfsTrace = handles.SniffTrace.MFS;
timeTrace = handles.SniffTrace.Timestamps;
for n = 1:numel(sniffs_of_interest)
    whichsniff = sniffs_of_interest(n);
    if isnan(handles.SniffsMFS(whichsniff,4))
        traceidx = handles.SniffsMFS(whichsniff,8);
        if abs(timeTrace(traceidx) - handles.SniffsMFS(whichsniff,1)) > 0.0025
            keyboard;
        end
        if mfsTrace(traceidx) <= 0
            f = find(mfsTrace(1:traceidx)>0,1,'last');
        else
            f = find(mfsTrace(traceidx:end)<0,1,'first') - 1 + traceidx;
        end
        if abs(traceidx - f) > 25
            set(handles.pointercarat,'XData',timeTrace(traceidx),'YData',mfsTrace(traceidx)+0.025);
            answer = questdlg('Keep the current detection?', ...
                'Adjust inhalation onset', ...
                'Yes','Pick','Delete','Yes');

            % Handle response
            switch answer
                case 'Yes'
                    handles.SniffsMFS(whichsniff,4) = timeTrace(traceidx);
                case 'Pick'
                    [x] = ginput(1); % select a pair of peak and valley
                    [~,idx] = min(abs(handles.SniffTrace.Timestamps - x(1)));
                    handles.SniffsMFS(whichsniff,4) = handles.SniffTrace.Timestamps(idx);
                case 'Delete'
                    handles.SniffsMFS(whichsniff-1,3) = handles.SniffsMFS(whichsniff+1,1);
                    handles.SniffsMFS(whichsniff,[4 7]) = -1;
%                 case 'Other'
%                     keyboard;
            end

        else
            idx = sort([f traceidx]);

            x = mfsTrace(idx(1):idx(2));
            y = timeTrace(idx(1):idx(2));
            xq = 0;
            % if there are any duplicates in x, interpolation will not work
            if numel(unique(x)) < numel(x)
                [~, tokeep] = unique(x,'stable');
                x = x(tokeep);
                y = y(tokeep);
            end
            handles.SniffsMFS(whichsniff,4) = interp1(x,y,xq);
        end
    end
    % keep track of the inhalation detected on the thermistor side as well
    if ~isnan(handles.SniffsMFS(whichsniff,4))
        thermInh = find( (handles.SniffsTS(:,1)>=handles.SniffsMFS(whichsniff,4)) & (handles.SniffsTS(:,1)<=handles.SniffsMFS(whichsniff,2)) );
        if isempty(thermInh)
            % case 1 - thermistor inhalation was before the MFS inhalation
            prevsniff = find(handles.SniffsMFS(1:whichsniff-1,4)>=0,1,'last');
            if ~isnan(handles.SniffsMFS(prevsniff,5))
                tmax = max(handles.SniffsMFS(whichsniff,4),handles.SniffsMFS(whichsniff,1));
                putative = find( (handles.SniffsTS(:,1)>handles.SniffsMFS(prevsniff,5)) & (handles.SniffsTS(:,1)<tmax) );
                if ~isempty(putative)
                    thermInh = putative;
                else
                    set(handles.pointercarat,'XData',handles.SniffsMFS(whichsniff,4),'YData',0.025);
                    % draw a line on the thernmistor axes as well
                    ylims = get(handles.SniffingFiltered,'YLim');
                    %set(handles.pointerline,'XData',handles.SniffsMFS(whichsniff,4)*[1 1],'YData',ylims);

                    answer = questdlg('Ignore thermistor inhalation detection?', ...
                        'Corresponding Thermistor peak', ...
                        'Yes','Pick','Delete','Yes');

                    % Handle response
                    switch answer
                        case 'Yes'
                            thermInh = nan;
                        case 'Pick'
                            [x] = ginput(1); % select a pair of peak and valley
                            [~,idx] = min(abs(handles.SniffTrace.Timestamps - x(1)));
                            handles.SniffsMFS(whichsniff,5) = handles.SniffTrace.Timestamps(idx);
                            handles.SniffsMFS(whichsniff,6) = mfsTrace(idx);
                            handles.SniffsMFS(whichsniff,7) = -1; % extra sniff detection that was missed in the thermistor
                            thermInh = nan;
                        case 'Delete'
                            handles.SniffsMFS(prevsniff,3) = handles.SniffsMFS(whichsniff+1,1);
                            handles.SniffsMFS(whichsniff,[4 7]) = -1;
                            thermInh = nan;
                        case 'Quit'
                            keyboard;
                    end
                end
            end
        end
        if isempty(thermInh)
            set(handles.pointercarat,'XData',handles.SniffsMFS(whichsniff,4),'YData',0.025);
            answer = questdlg('Ignore thermistor inhalation detection?', ...
                'Corresponding Thermistor peak', ...
                'Yes','Pick','Delete','Yes');

            % Handle response
            switch answer
                case 'Yes'
                    thermInh = nan;
                case 'Pick'
                    [x] = ginput(1); % select a pair of peak and valley
                    [~,idx] = min(abs(handles.SniffTrace.Timestamps - x(1)));
                    handles.SniffsMFS(whichsniff,5) = handles.SniffTrace.Timestamps(idx);
                    handles.SniffsMFS(whichsniff,6) = mfsTrace(idx);
                    handles.SniffsMFS(whichsniff,7) = -1; % extra sniff detection that was missed in the thermistor
                    thermInh = nan;
                case 'Delete'
                    handles.SniffsMFS(prevsniff,3) = handles.SniffsMFS(whichsniff+1,1);
                    handles.SniffsMFS(whichsniff,[4 7]) = -1;
                    thermInh = nan;
                case 'Quit'
                    keyboard;
            end
        end
        if ~isnan(thermInh(1))
            handles.SniffsMFS(whichsniff,5) = handles.SniffsTS(thermInh(1),1);
            handles.SniffsMFS(whichsniff,6) = mfsTrace(handles.SniffsTS(thermInh(1),8));
        elseif handles.SniffsMFS(whichsniff,7) ~= -1
            % handles.SniffsMFS(whichsniff,5:6) = nan; 
        end
    end
end
% delete flagged sniffs
handles.SniffsMFS(find(handles.SniffsMFS(:,4)==-1),:) = [];
set(handles.pointercarat,'XData',[],'YData',[]);
UpdatePeakValleyPlots(hObject, eventdata, handles);
if ~isfield(handles,'roiPosition')
    delete(roi);
end
guidata(hObject, handles);

function DeleteZeroCrossings_Callback(hObject, eventdata, handles)
axes(handles.MFS);
if ~isfield(handles,'roiPosition')
    roi = drawrectangle;
    sniffs_of_interest = intersect(...
                        find(handles.SniffsMFS(:,1) >= roi.Position(1)), ...
                            find(handles.SniffsMFS(:,2) <= (roi.Position(1) + roi.Position(3)) ) );
else
    sniffs_of_interest = intersect(...
                        find(handles.SniffsMFS(:,1) >= handles.roiPosition(1)), ...
                            find(handles.SniffsMFS(:,2) <= (handles.roiPosition(2)) ) );
end

handles.SniffsMFS(sniffs_of_interest,4:6) = nan;
UpdatePeakValleyPlots(hObject, eventdata, handles);
if ~isfield(handles,'roiPosition')
    delete(roi);
end
guidata(hObject, handles);


% --- Executes on button press in RedoStretch.
function RedoStretch_Callback(hObject, eventdata, handles)
SniffTag = handles.primarySniffTS;
AxesTag = handles.primaryAxes; %'MFS2Therm'; % 'SniffingFiltered'
TraceTag = handles.primaryRawSniffTrace; %'mfs2thermTrace'; % 'rawTrace'
SniffTraceTag = handles.primarySniffTrace; %'MFS2Therm'; %'Filtered'

axes(handles.(AxesTag));
roi = drawrectangle;
indices_of_interest = intersect(...
                    find(handles.SniffTrace.Timestamps >= roi.Position(1)), ...
                        find(handles.SniffTrace.Timestamps <= (roi.Position(1) + roi.Position(3)) ) );

redo = 1;
while redo
    redo = redo + 1;
    sdnew = redo*str2double(handles.SDfactor.String);

    stretchTimeStamps = handles.(TraceTag).XData(indices_of_interest)';
    stretchThermistor = handles.(TraceTag).YData(indices_of_interest)';
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
            'YData',handles.SniffTrace.(SniffTraceTag)(thisStretchSniffs(:,8)) );
        set(handles.valleysNew, ...
            'XData',handles.SniffTrace.Timestamps(thisStretchSniffs(:,9)),...
            'YData',handles.SniffTrace.(SniffTraceTag)(thisStretchSniffs(:,9)) );

        % keep the new detections
        answer = questdlg('Keep the new detections?', ...
            'Redo stretch', ...
            'Yes','Redo','Cancel','Yes');

        % Handle response
        switch answer
            case 'Yes'
                % integrate into the current set of sniffs
                whichpeaks = find( (handles.(SniffTag)(:,8)>=indices_of_interest(1)) & (handles.(SniffTag)(:,8)<=indices_of_interest(end)) );
                whichvalls = find( (handles.(SniffTag)(:,9)>=indices_of_interest(1)) & (handles.(SniffTag)(:,9)<=indices_of_interest(end)) );
                rows2replace = unique(vertcat(whichpeaks,whichvalls));

                if isempty(rows2replace)
                    % integrate in between sniffs
                    prevsniff = find(handles.(SniffTag)(:,1)<thisStretchSniffs(1,1),1,'last');
                    nextsniff = find(handles.(SniffTag)(:,1)>thisStretchSniffs(end,1),1,'first');
                    if isempty(prevsniff)
                        chunkC = handles.(SniffTag)(nextsniff:end,:);
                        chunkB = thisStretchSniffs;
                        if strcmp(handles.primarySniffTS,'SniffsMFS')
                            chunkB(:,4:7) = nan;
                        else
                            chunkB(:,4:7) = repmat(handles.SniffsTS(rows2replace(1),4:7),size(thisStretchSniffs,1),1);
                        end
                        chunkB(1,3) = chunkC(1,1);
                        handles.(SniffTag) = [chunkB; chunkC];
                    elseif  nextsniff - prevsniff == 1 && size(thisStretchSniffs,1) == 1
                        chunkA = handles.(SniffTag)(1:prevsniff,:);
                        chunkA(end,3) = thisStretchSniffs(1,1);
                        chunkC = handles.(SniffTag)(nextsniff:end,:);
                        chunkB = thisStretchSniffs;
                        if strcmp(handles.primarySniffTS,'SniffsMFS')
                            chunkB(:,4:7) = nan;
                        else
                            chunkB(:,4:7) = repmat(handles.SniffsTS(rows2replace(1),4:7),size(thisStretchSniffs,1),1);
                        end
                        chunkB(1,3) = chunkC(1,1);
                        handles.(SniffTag) = [chunkA; chunkB; chunkC];
                    end
                    redo = 0;
                else
                    if numel(rows2replace) > 1
                        % is the last detection before a whole sniff?
                        if thisStretchSniffs(end,2) < handles.(SniffTag)(rows2replace(end),1)
                            rows2replace(end,:) = [];
                        end
                    end

                    chunkA = handles.(SniffTag)(1:(rows2replace(1)-1),:);
                    chunkC = handles.(SniffTag)((rows2replace(end)+1):end,:);
                    chunkB = thisStretchSniffs;
                    if strcmp(handles.primarySniffTS,'SniffsMFS')
                        chunkB(:,4:7) = nan;
                    else
                        chunkB(:,4:7) = repmat(handles.SniffsTS(rows2replace(1),4:7),size(thisStretchSniffs,1),1);
                    end
                    chunkB(1,[1 8]) = handles.(SniffTag)(rows2replace(1),[1 8]);
                    chunkB(end,3) = chunkC(1,1);
                    handles.(SniffTag) = [chunkA; chunkB; chunkC];

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
% --- Executes on button press in FlagStretch.
function FlagStretch_Callback(hObject, eventdata, handles)
roi = drawrectangle;

SniffTag = handles.primarySniffTS;

sniffs_to_flag = intersect(...
                    find(handles.(SniffTag)(:,1) >= roi.Position(1)), ...
                        find(handles.(SniffTag)(:,1) <= (roi.Position(1) + roi.Position(3)) ) );
handles.(SniffTag)(sniffs_to_flag,8:9) = -abs(handles.(SniffTag)(sniffs_to_flag,8:9));

UpdatePeakValleyPlots(hObject, eventdata, handles);
                    
% axes(handles.SniffingRaw);
% set(handles.rawTrace,'XData', handles.SniffTrace.Timestamps, 'YData', handles.SniffTrace.Raw);
delete(roi);
guidata(hObject, handles);

% ===================================================================================
% --- Executes on button press in AddSniff.
function SplitSniff_Callback(hObject, eventdata, handles)
% hObject    handle to AddSniff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.(handles.primaryAxes));
SniffTag = handles.primarySniffTS;
zoom off
[x,y] = ginput(2); % select a pair of peak and valley

[~,idx1] = min(abs(handles.SniffTrace.Timestamps - x(1)));
[~,idx2] = min(abs(handles.SniffTrace.Timestamps - x(2)));
handles.newSniff = [];
handles.newSniff(1,[2 9]) = [handles.SniffTrace.Timestamps(idx1) idx1];
handles.newSniff(2,[1 8]) = [handles.SniffTrace.Timestamps(idx2) idx2];

% integrate into the current set of sniffs
whichsniff = find(handles.(SniffTag)(:,1)<=handles.newSniff(1,2),1,'last');
handles.newSniff(2,[2 3:7 9]) = handles.(SniffTag)(whichsniff,[2 3:7 9]);
handles.(SniffTag)(whichsniff,[2 9]) = handles.newSniff(1,[2 9]);
handles.(SniffTag)(whichsniff,3) = handles.newSniff(2,1);
handles.newSniff(1,:) = [];
handles.(SniffTag) = [handles.(SniffTag)(1:whichsniff,:); handles.newSniff; handles.(SniffTag)((whichsniff+1):end,:)];
guidata(hObject, handles);
UpdatePeakValleyPlots(hObject, eventdata, handles); 

% --- Executes on button press in AddSniff.
function AddSniff_Callback(hObject, eventdata, handles)
% hObject    handle to AddSniff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.(handles.primaryAxes));
SniffTag = handles.primarySniffTS;
zoom off
[x,y] = ginput(2); % select a pair of peak and valley

[~,idx1] = min(abs(handles.SniffTrace.Timestamps - x(1)));
[~,idx2] = min(abs(handles.SniffTrace.Timestamps - x(2)));
handles.newSniff(1,[1 8]) = [handles.SniffTrace.Timestamps(idx1) idx1];
handles.newSniff(1,[2 9]) = [handles.SniffTrace.Timestamps(idx2) idx2];

% integrate into the current set of sniffs
whichsniff = find(handles.(SniffTag)(:,1)<=handles.newSniff(1,1),1,'last');
handles.newSniff(1,3:7) = handles.(SniffTag)(whichsniff,3:7);
handles.(SniffTag)(whichsniff,3) = handles.newSniff(1,1);
handles.(SniffTag) = [handles.(SniffTag)(1:whichsniff,:); handles.newSniff; handles.(SniffTag)((whichsniff+1):end,:)];
guidata(hObject, handles);
UpdatePeakValleyPlots(hObject, eventdata, handles); 

% --- Executes on button press in TempSave.
function TempSave_Callback(hObject, eventdata, handles)
if strcmp(handles.primarySniffTS,'SniffsMFS')
    MFSDetectionThreshold = str2double(handles.SDfactor.String);
    CuratedMFS_SniffTimestamps = handles.SniffsMFS;
    lastMFSTimestamp = floor(handles.SniffingFiltered.XLim(2));
    save(handles.WhereSession.String,'MFSDetectionThreshold','CuratedMFS_SniffTimestamps','lastMFSTimestamp','-append');
    disp(['saved MFS sniffs! @ ',num2str(lastMFSTimestamp), ' seconds']);
else
    SniffDetectionThreshold = str2double(handles.SDfactor.String);
    Curated_SniffTimestamps = handles.SniffsTS;
    lastTimestamp = floor(handles.SniffingFiltered.XLim(2));
    save(handles.WhereSession.String,'SniffDetectionThreshold','Curated_SniffTimestamps','lastTimestamp','-append');
    disp(['saved Thermistor sniffs! @ ',num2str(lastTimestamp), ' seconds']);
end

% --- Executes on button press in SaveSniffs.
function SaveSniffs_Callback(hObject, eventdata, handles)
if strcmp(handles.primarySniffTS,'SniffsMFS')
    MFSDetectionThreshold = str2double(handles.SDfactor.String);
    CuratedMFS_SniffTimestamps = handles.SniffsMFS;
    CuratedMFSSniffTimestamps = handles.SniffsMFS;
    lastMFSTimestamp = floor(handles.SniffingFiltered.XLim(2));
    save(handles.WhereSession.String,'MFSDetectionThreshold','CuratedMFS_SniffTimestamps','CuratedMFSSniffTimestamps','lastMFSTimestamp','-append');
    disp(['saved sniffs! @ ',num2str(lastMFSTimestamp), ' seconds']);
else
    SniffDetectionThreshold = str2double(handles.SDfactor.String);
    Curated_SniffTimestamps = handles.SniffsTS;
    CuratedSniffTimestamps = handles.SniffsTS;
    lastTimestamp = floor(handles.SniffingFiltered.XLim(2));
    save(handles.WhereSession.String,'SniffDetectionThreshold','Curated_SniffTimestamps','CuratedSniffTimestamps','lastTimestamp','-append');
    disp(['saved sniffs! @ ',num2str(lastTimestamp), ' seconds']);
end

% --- Executes on button press in RecoverTemp.
function RecoverTemp_Callback(hObject, eventdata, handles)
if strcmp(handles.primarySniffTS,'SniffsMFS')
    load(handles.WhereSession.String,'CuratedMFS_SniffTimestamps','lastMFSTimestamp');
    if ~exist('lastMFSTimestamp','var')
        TSdone = input('Timestamp until which data was curated?');
    else
        TSdone = lastMFSTimestamp;
    end
    whichidx1 = find(CuratedMFS_SniffTimestamps(:,1)>=TSdone,1,'first');
    whichidx2 = find(handles.SniffsMFS(:,1)>=TSdone,1,'first');
    if isempty(whichidx1)
        handles.SniffsMFS = CuratedMFS_SniffTimestamps;
    elseif isequal(CuratedMFS_SniffTimestamps(whichidx1,1:3),handles.SniffsMFS(whichidx2,1:3))
        handles.SniffsMFS = [CuratedMFS_SniffTimestamps(1:(whichidx1-1),:); handles.SniffsMFS(whichidx2:end,:)];
    else
        keyboard;
    end
else
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
end
guidata(hObject, handles); 
UpdatePeakValleyPlots(hObject, eventdata, handles);

%% keyboard shortcuts
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
    case {'uparrow'}
        handles.roiPosition = get(handles.SniffingFiltered,'XLim');
        GetZeroCrossings_Callback(hObject, eventdata, handles);
    case {'downarrow'}
        handles.roiPosition = get(handles.SniffingFiltered,'XLim');
        DeleteZeroCrossings_Callback(hObject, eventdata, handles);
    case 'rightarrow'
        NextStretch_Callback(hObject, eventdata, handles);
    case 'leftarrow'
        PreviousStretch_Callback(hObject, eventdata, handles);
    case 'm'
        axes(handles.(handles.primaryAxes));
        ModifyPoint_Callback(hObject, eventdata, handles);
    case 'n'
        axes(handles.(handles.primaryAxes));
        AddSniff_Callback(hObject, eventdata, handles);
    case 'b'
        axes(handles.(handles.primaryAxes));
        SplitSniff_Callback(hObject, eventdata, handles);
    case 'r'
        axes(handles.(handles.primaryAxes));
        RedoStretch_Callback(hObject, eventdata, handles);
    case 'o'
        % Call a function or perform an action for 'o' key
        disp('Open action triggered!');
        % call open_file_function(handles);
        % Add more cases for other keys
%     case 'z'
%         axes(handles.SniffingFiltered);
%         zoom yon
%         guidata(hObject, handles);

end
