function varargout = ProcessSniffTimeStamps_GUI(varargin)
% PROCESSSNIFFTIMESTAMPS_GUI MATLAB code for ProcessSniffTimeStamps_GUI.fig
%      PROCESSSNIFFTIMESTAMPS_GUI, by itself, creates a new PROCESSSNIFFTIMESTAMPS_GUI or raises the existing
%      singleton*.
%
%      H = PROCESSSNIFFTIMESTAMPS_GUI returns the handle to a new PROCESSSNIFFTIMESTAMPS_GUI or the handle to
%      the existing singleton*.
%
%      PROCESSSNIFFTIMESTAMPS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCESSSNIFFTIMESTAMPS_GUI.M with the given input arguments.
%
%      PROCESSSNIFFTIMESTAMPS_GUI('Property','Value',...) creates a new PROCESSSNIFFTIMESTAMPS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProcessSniffTimeStamps_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProcessSniffTimeStamps_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProcessSniffTimeStamps_GUI

% Last Modified by GUIDE v2.5 10-Feb-2025 09:02:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcessSniffTimeStamps_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcessSniffTimeStamps_GUI_OutputFcn, ...
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


% --- Executes just before ProcessSniffTimeStamps_GUI is made visible.
function ProcessSniffTimeStamps_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProcessSniffTimeStamps_GUI (see VARARGIN)

% Choose default command line output for ProcessSniffTimeStamps_GUI
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
handles.rawTrace    = plot(nan,nan,'color',Plot_Colors('b'));
handles.peaksRaw    = plot(nan,nan,'ok','MarkerSize',4); 
handles.valleysRaw  = plot(nan,nan,'or','MarkerSize',4);

axes(handles.SniffingFiltered);
hold on;
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

%if exist(handles.WhereSession.String)==2
    LoadSession_Callback(hObject, eventdata, handles);
%end

% UIWAIT makes ProcessSniffTimeStamps_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ProcessSniffTimeStamps_GUI_OutputFcn(hObject, eventdata, handles) 
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
        load(handles.WhereSession.String, 'AllSniffs', 'RespirationData'); % AllSniffs: nx13, RespirationData: t x 3(ts, raw, filt)
        handles.SessionLength = RespirationData(end,1);

        try
            load(fullfile(fileparts(handles.WhereSession.String), 'quickprocessOdorTTLs.mat'));
        catch
            warning('no TTLs found');
        end

        % make equivalent long traces for plotting
        handles.SniffTrace.Timestamps   = RespirationData(:,1);
        handles.OdorLocationTrace       = [];
        handles.SniffTrace.Raw          = RespirationData(:,2); % unfiltered thermistor trace
        handles.SniffTrace.Filtered     = RespirationData(:,2); % filtered thermistor trace
        
        handles.SniffsTS(:,1:2) = AllSniffs(:,1:2);
        handles.SniffsTS(:,8:9) = AllSniffs(:,11:12);
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

% plot both filtered and raw traces, and detcted timestamps
axes(handles.SniffingFiltered);
set(handles.filtTrace,'XData', handles.SniffTrace.Timestamps, 'YData', handles.SniffTrace.Filtered);
zoom off
% overlay detected timestamps on the plot
set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,8)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTS(:,8)));
set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,9)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTS(:,9)));

set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)], ...
    'YTick', []);
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

axes(handles.SniffingRaw);
set(handles.rawTrace,'XData', handles.SniffTrace.Timestamps, 'YData', handles.SniffTrace.Raw);
zoom off
% overlay detected timestamps on the plot
set(handles.peaksRaw,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,8)),...
    'YData',handles.SniffTrace.Raw(handles.SniffsTS(:,8)));
set(handles.valleysRaw,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(:,9)),...
    'YData',handles.SniffTrace.Raw(handles.SniffsTS(:,9)));
set(gca,'XLim',handles.SniffTrace.Timestamps(1)+[0 str2double(handles.WindowSize.String)], ...
     'YTick', [],  'XTick', []);

if strcmp(handles.datamode,'cid')
    line(handles.SniffTrace.Timestamps([1 end]),[0 0],'color','k','Linestyle',':');
end

axes(handles.SniffingFiltered);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in RefindPeaks.
function RefindPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to RefindPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sdnew = 2*str2double(handles.SDfactor.String);

switch handles.datamode
    case {'cid', 'onlyEphys'}
        % re-detect peaks and valleys with lower peak prominence
        [handles.SniffsTSnew, sniffindices] = ProcessThermistorData([handles.SniffTrace.Timestamps handles.SniffTrace.Filtered],'SDfactor',sdnew);
        handles.SniffsTSnew(:,8:9) = sniffindices(:,1:2);
        handles.SniffsTSnew(find(isnan(handles.SniffsTSnew(:,8))),:) = [];
        handles.SniffsTSnew(find(handles.SniffsTSnew(:,8)<=0),:) = [];

        % find trace indices that correspond to detected sniff timestamps
        for n = 1:size(handles.SniffsTSnew,1) % every new sniff detected

            % was this already detected
            f = find(abs(handles.SniffsTS(:,1)-handles.SniffsTSnew(n,1))<0.004,1,'first');
            if ~isempty(f) & abs(handles.SniffsTS(f,2)-handles.SniffsTSnew(n,2))<0.004
                handles.SniffsTSnew(n,8:9) = NaN;
            end
        end
    case 'smellocator'
        % re-detect peaks and valleys with lower peak prominence
        load(handles.WhereSession.String,'Traces','TrialInfo');
        if ~isfield(Traces,'OdorLocation')
            Traces.OdorLocation     = Traces.Motor;
        end
        [SniffTimeStamps] = ...
            TrialWiseSniffs(TrialInfo,Traces,'SDfactor',sdnew); % [sniffstart sniffstop nextsniff odorlocation sniffslope stimstate trialID]
        % remove overlapping sniffs
        handles.SniffsTSnew = SniffTimeStamps(find(SniffTimeStamps(:,end)>0),:);

        % find trace indices that correspond to detected sniff timestamps
        for n = 1:size(handles.SniffsTSnew,1)

            % was this already detected
            f = find(abs(handles.SniffsTS(:,1)-handles.SniffsTSnew(n,1))<0.004,1,'first');
            if ~isempty(f) & abs(handles.SniffsTS(f,2)-handles.SniffsTSnew(n,2))<0.004
                handles.SniffsTSnew(n,8:9) = NaN;
            else
                % inhalation start
                [~,idx] = min(abs(handles.SniffTrace.Timestamps - handles.SniffsTSnew(n,1)));
                if abs(handles.SniffTrace.Timestamps(idx) - handles.SniffsTSnew(n,1)) < 0.004
                    handles.SniffsTSnew(n,8) = idx;
                end
                % inhalation end
                [~,idx] = min(abs(handles.SniffTrace.Timestamps - handles.SniffsTSnew(n,2)));
                if abs(handles.SniffTrace.Timestamps(idx) - handles.SniffsTSnew(n,2)) < 0.004
                    handles.SniffsTSnew(n,9) = idx;
                end
            end
        end
end



% collate the list and make sure peaks fall in order
[handles] = collatesniffs(handles,0);
% update plots
UpdatePeakValleyPlots(hObject, eventdata, handles);
handles.KeepNewSD.BackgroundColor = [0 0.94 0];
handles.KeepOldSD.BackgroundColor = [0 0.94 0];
guidata(hObject, handles);
%uiwait;
%guidata(hObject, handles);


function [handles] = collatesniffs(handles,whichmode)

% SniffTSnew has new peaks/valleys detected with the new lower user
% selected threshold
% many of these peaks/valleys are redundant with those detected with the
% standard threshold (2.5 SD)
% redundant new detections are flagged by swapping col 8,9 values with nan
if isfield(handles,'SniffsTSnew')
    newsniffs = handles.SniffsTSnew(~isnan(handles.SniffsTSnew(:,8)),:);
    % swap col 8,9 (trace indices) to -ve to know in the pooled set which
    % sniffs came for new detections
    newsniffs(:,8:9) = -newsniffs(:,8:9);
else
    newsniffs = [];
end
% pool old sniffs and new detections
pooledsniffs = vertcat(handles.SniffsTS,newsniffs);
pooledsniffs = sortrows(pooledsniffs,1);

% for every new detection - reconcile with previous detections
while any(pooledsniffs(:,9)<0)
    f = find(pooledsniffs(:,9)<0,1,'first');
    if numel(f)<2
        break;
    end
%     newLims = pooledsniffs(f,1) + [-floor(str2double(handles.WindowSize.String)/2) ceil(str2double(handles.WindowSize.String)/2)];
%     set(handles.SniffingRaw,'XLim',newLims);
%     set(handles.SniffingFiltered,'XLim',newLims);
%     handles.Scroller.Value = pooledsniffs(f,1)/(handles.SessionLength - str2double(handles.WindowSize.String));
    
    if pooledsniffs(f,1)>pooledsniffs(f-1,2) & pooledsniffs(f,1)<pooledsniffs(f-1,3)
        % new sniff start and end falls within a previous sniff 
        k = find(abs(pooledsniffs(f:end,3)-pooledsniffs(f-1,3))<0.01); % previous sniffs end must match another sniff in the new pooled set
        if ~isempty(k)
            k = k + f - 1;
            pooledsniffs(f-1,3) = pooledsniffs(f,1); % change the sniff end for the previous sniff
            pooledsniffs(f:k,9) = -pooledsniffs(f:k,9); % un-negate col 9 for new detections that can be reconciled 
        else
            keyboard;
        end
    elseif pooledsniffs(f,1)>pooledsniffs(f-1,2) & pooledsniffs(f,1)==pooledsniffs(f-1,3)
        % if the new detected sniff overlaps perfect with previous sniff end 
        % pooledsniffs(f-1,3) = pooledsniffs(f,1); % superfluous
        pooledsniffs(f,9) = -pooledsniffs(f,9);        
        k = find(abs(pooledsniffs(f+1:end,1)-pooledsniffs(f-1,3))<0.01); % previous sniffs end must match another sniff in the non-pooled set
        if ~isempty(k)
            k = k + f + 1 - 1;
            keyboard;
            pooledsniffs(k,8:9) = -abs(pooledsniffs(k,8:9)); 
        else
            %keyboard;
        end
    elseif abs(pooledsniffs(f,1)-pooledsniffs(f-1,1))<0.01
        % just high overlap but not an exact match
        k = find(abs(pooledsniffs(f:end,3)-pooledsniffs(f-1,3))<0.01); % find matching ends
        if ~isempty(k)
            k = k + f - 1;
            m = find(handles.SniffsTS(:,1)==pooledsniffs(f-1,1));
            pooledsniffs(f:k,9) = -pooledsniffs(f:k,9);
            pooledsniffs(f-1,:) = [];
            handles.SniffsTS(m,8) = -abs(handles.SniffsTS(m,8));
            
%             % overlay detected timestamps on the plot
%             set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),8)),...
%             'YData',handles.SniffTrace.Filtered(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),8)));
%             set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),9)),...
%             'YData',handles.SniffTrace.Filtered(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),9)));
            
            
        else
            keyboard;
        end
    else
        keyboard;
    end
end

if whichmode
    % for saving
    % sanity check 1: any duplicates?
    while any(abs(diff(pooledsniffs(:,1)))<0.01)
        f = find(abs(diff(pooledsniffs(:,1)))<0.01,1,'first');
        if pooledsniffs(f,1) == pooledsniffs(f+1,1) & pooledsniffs(f,2) == pooledsniffs(f+2,2)
            pooledsniffs(f,:) = [];
        elseif pooledsniffs(f,1) == pooledsniffs(f+1,1) & pooledsniffs(f,2) == pooledsniffs(f+3,2)
            pooledsniffs(f,:) = [];
        else
            keyboard;
        end
    end
    
    % sanity check 2: any gaps?
    ignoreme = 0;
    while any(abs(pooledsniffs(1:end-1,3)-pooledsniffs(2:end,1))>0.01) & ~ignoreme
        f = find(abs(pooledsniffs(1:end-1,3)-pooledsniffs(2:end,1))>0.01,1,'first');
        if pooledsniffs(f+1,1)<pooledsniffs(f,3) & pooledsniffs(f+1,1)>pooledsniffs(f,2) & pooledsniffs(f,3)==pooledsniffs(f+1,3)
            pooledsniffs(f,3) = pooledsniffs(f+1,1);
        elseif abs(pooledsniffs(f+1,2)-pooledsniffs(f,2))<=0.005 & abs(pooledsniffs(f+1,3)-pooledsniffs(f,3))<=0.005
            keyboard;
            %todelete = vertcat(todelete,f(x)+1);
            %pooledsniffs(f+1,:) = [];
        elseif pooledsniffs(f+[1:2],8) < 0
            pooledsniffs(f+[1:2],:) = [];
        else
            keyboard;
            % pooledsniffs(f,3) = pooledsniffs(f+1,1);
        end
    end

    % remove negative indices
    pooledsniffs(:,8:9) = abs(pooledsniffs(:,8:9));

    % sanity check 3: indices match timestamps?
    if any(pooledsniffs(:,1)~=handles.SniffTrace.Timestamps(pooledsniffs(:,8)))
        if ~any(abs(pooledsniffs(:,1)-handles.SniffTrace.Timestamps(pooledsniffs(:,8)))>=0.002)
            f = find((pooledsniffs(:,1)~=handles.SniffTrace.Timestamps(pooledsniffs(:,8))));
            for x = 1:numel(f)
                if pooledsniffs(f(x)-1,3) == pooledsniffs(f(x),1)
                    pooledsniffs(f(x)-1,3) = handles.SniffTrace.Timestamps(pooledsniffs(f(x),8));
                else
                    keyboard;
                end
                pooledsniffs(f(x),1) = handles.SniffTrace.Timestamps(pooledsniffs(f(x),8));
            end
        else
            keyboard;
        end
    end
    if any(pooledsniffs(:,2)~=handles.SniffTrace.Timestamps(pooledsniffs(:,9)))
        if ~any(abs(pooledsniffs(:,2)-handles.SniffTrace.Timestamps(pooledsniffs(:,9)))>=0.002)
            f = find((pooledsniffs(:,2)~=handles.SniffTrace.Timestamps(pooledsniffs(:,9))));
            pooledsniffs(f,2) = handles.SniffTrace.Timestamps(pooledsniffs(f,9));
        else
            keyboard;
        end
    end
    
    % flag any sniffs in the flagged strectch
    bad_indices = find(isnan(handles.rawTrace.YData));
    if ~isempty(bad_indices)
        segment_idx(:,1) = [bad_indices(1) bad_indices(find(diff(bad_indices)~=1)+1)];
        segment_idx(:,2) = [bad_indices(find(diff(bad_indices)~=1)) bad_indices(end)];
%         terribleTimestamps = [];
%         segments = numel(find(diff(bad_indices)>1)) + 1;
%         if segments == 1
%             terribleTimestamps(1,:) = handles.rawTrace.XData(bad_indices([1 end]));
%             % flag any sniffs that are in this period
%             iffysniffs = find(pooledsniffs(:,1)>terribleTimestamps(1,1),1,'first') : ...
%                 find(pooledsniffs(:,1)<terribleTimestamps(1,2),1,'last');
%             pooledsniffs(iffysniffs,10) = -1;
%         else
%             
%         end
        for segment = 1:size(segment_idx,1)
            terribleTimestamps = [];
            terribleTimestamps(1,:) = handles.rawTrace.XData(segment_idx(segment,:));
            % flag any sniffs that are in this period
            iffysniffs = find(pooledsniffs(:,1)>terribleTimestamps(1,1),1,'first') : ...
                find(pooledsniffs(:,1)<terribleTimestamps(1,2),1,'last');
            pooledsniffs(iffysniffs,10) = -1;
        end
    end
        
    % update plot
    set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(pooledsniffs(:,8)),...
        'YData',handles.SniffTrace.Filtered(pooledsniffs(:,8)));
    set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(pooledsniffs(:,9)),...
        'YData',handles.SniffTrace.Filtered(pooledsniffs(:,9)));

    set(handles.peaksNew, 'XData',[], 'YData',[]);
    set(handles.valleysNew,  'XData',[], 'YData',[]);

    % append timestamps to processed file
    if strcmp(handles.datamode,'smellocator_raw')
        PassiveSniffDetectionThreshold = str2double(handles.SDfactor.String);
        CuratedPassiveSniffTimestamps = pooledsniffs;

        save(handles.ProcessedSession,'PassiveSniffDetectionThreshold','CuratedPassiveSniffTimestamps','-append');
    else
        SniffDetectionThreshold = str2double(handles.SDfactor.String);
        CuratedSniffTimestamps = pooledsniffs;

        save(handles.WhereSession.String,'SniffDetectionThreshold','CuratedSniffTimestamps','-append');
    end
    disp('saved sniffs!')
end

function [handles] = collatesniffs_old(handles,whichmode)

% SniffTSnew has new peaks/valleys detected with the new lower user
% selected threshold
% many of these peaks/valleys are redundant with those detected with the
% standard threshold (2.5 SD)
% redundant new detections are flagged by swapping col 8,9 values with nan
newsniffs = handles.SniffsTSnew(~isnan(handles.SniffsTSnew(:,8)),:);
% swap col 8,9 (trace indices) to -ve to know in the pooled set which
% sniffs came for new detections
newsniffs(:,8:9) = -newsniffs(:,8:9);
% pool old sniffs and new detections
pooledsniffs = vertcat(handles.SniffsTS,newsniffs);
pooledsniffs = sortrows(pooledsniffs,1);

% for every new detection - reconcile with previous detections
while any(pooledsniffs(:,9)<0)
    f = find(pooledsniffs(:,9)<0,1,'first');
    if numel(f)<2
        break;
    end
%     newLims = pooledsniffs(f,1) + [-floor(str2double(handles.WindowSize.String)/2) ceil(str2double(handles.WindowSize.String)/2)];
%     set(handles.SniffingRaw,'XLim',newLims);
%     set(handles.SniffingFiltered,'XLim',newLims);
%     handles.Scroller.Value = pooledsniffs(f,1)/(handles.SessionLength - str2double(handles.WindowSize.String));
    
    if pooledsniffs(f,1)>pooledsniffs(f-1,2) & pooledsniffs(f,1)<pooledsniffs(f-1,3)
        % new sniff start and end falls within a previous sniff 
        k = find(abs(pooledsniffs(f:end,3)-pooledsniffs(f-1,3))<0.01); % previous sniffs end must match another sniff in the new pooled set
        if ~isempty(k)
            k = k + f - 1;
            pooledsniffs(f-1,3) = pooledsniffs(f,1); % change the sniff end for the previous sniff
            pooledsniffs(f:k,9) = -pooledsniffs(f:k,9); % un-negate col 9 for new detections that can be reconciled 
        else
            keyboard;
        end
    elseif pooledsniffs(f,1)>pooledsniffs(f-1,2) & pooledsniffs(f,1)==pooledsniffs(f-1,3)
        % if the new detected sniff overlaps perfect with previous sniff end 
        % pooledsniffs(f-1,3) = pooledsniffs(f,1); % superfluous
        pooledsniffs(f,9) = -pooledsniffs(f,9);        
        k = find(abs(pooledsniffs(f+1:end,1)-pooledsniffs(f-1,3))<0.01); % previous sniffs end must match another sniff in the non-pooled set
        if ~isempty(k)
            k = k + f + 1 - 1;
            keyboard;
            pooledsniffs(k,8:9) = -abs(pooledsniffs(k,8:9)); 
        else
            %keyboard;
        end
    elseif abs(pooledsniffs(f,1)-pooledsniffs(f-1,1))<0.01
        % just high overlap but not an exact match
        k = find(abs(pooledsniffs(f:end,3)-pooledsniffs(f-1,3))<0.01); % find matching ends
        if ~isempty(k)
            k = k + f - 1;
            m = find(handles.SniffsTS(:,1)==pooledsniffs(f-1,1));
            pooledsniffs(f:k,9) = -pooledsniffs(f:k,9);
            pooledsniffs(f-1,:) = [];
            handles.SniffsTS(m,8) = -abs(handles.SniffsTS(m,8));
            
%             % overlay detected timestamps on the plot
%             set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),8)),...
%             'YData',handles.SniffTrace.Filtered(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),8)));
%             set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),9)),...
%             'YData',handles.SniffTrace.Filtered(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),9)));
            
            
        else
            keyboard;
        end
    else
        keyboard;
    end
end

if whichmode
    % for saving
    % sanity check 1: any duplicates?
    if any(abs(diff(pooledsniffs(:,1)))<0.01)
        f = find(abs(diff(pooledsniffs(:,1)))<0.01);
        if pooledsniffs(f,1) == pooledsniffs(f+1,1) & pooledsniffs(f,2) == pooledsniffs(f+2,2)
            pooledsniffs(f,:) = [];
            if any(abs(diff(pooledsniffs(:,1)))<0.01)
                keyboard;
            end
        else
            keyboard;
        end
    end
    % sanity check 2: any gaps?
    if any(abs(pooledsniffs(1:end-1,3)-pooledsniffs(2:end,1))>0.01)
        todelete = [];
        f = find(abs(pooledsniffs(1:end-1,3)-pooledsniffs(2:end,1))>0.01);
        for x = 1:numel(f)
            if pooledsniffs(f(x)+1,1)<pooledsniffs(f(x),3) & pooledsniffs(f(x)+1,1)>pooledsniffs(f(x),2) & pooledsniffs(f(x),3)==pooledsniffs(f(x)+1,3)
                pooledsniffs(f(x),3) = pooledsniffs(f(x)+1,1);
            elseif abs(pooledsniffs(f(x)+1,2)-pooledsniffs(f(x),2))<=0.005 & abs(pooledsniffs(f(x)+1,3)-pooledsniffs(f(x),3))<=0.005
                keyboard;
                %todelete = vertcat(todelete,f(x)+1);
            end
        end
        if ~isempty(todelete)
            pooledsniffs(todelete,:) = [];
        end
        if any(abs(pooledsniffs(1:end-1,3)-pooledsniffs(2:end,1))>0.01)
            keyboard;
        end
    end

    % remove negative indices
    pooledsniffs(:,8:9) = abs(pooledsniffs(:,8:9));

    % sanity check 3: indices match timestamps?
    if any(pooledsniffs(:,1)~=handles.SniffTrace.Timestamps(pooledsniffs(:,8))) ...
            || ...
            any(pooledsniffs(:,2)~=handles.SniffTrace.Timestamps(pooledsniffs(:,9)))
        keyboard;
    end
    
    % flag any sniffs in the flagged strectch
    bad_indices = find(isnan(handles.rawTrace.YData));
    if ~isempty(bad_indices)
        terribleTimestamps = [];
        segments = numel(find(diff(bad_indices)>1)) + 1;
        if segments == 1
            terribleTimestamps(1,:) = handles.rawTrace.XData(bad_indices([1 end]));
            % flag any sniffs that are in this period
            iffysniffs = find(pooledsniffs(:,1)>terribleTimestamps(1,1),1,'first') : ...
                find(pooledsniffs(:,1)<terribleTimestamps(1,2),1,'last');
            pooledsniffs(iffysniffs,10) = -1;
        else
            
        end
    end
        
    % update plot
    set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(pooledsniffs(:,8)),...
        'YData',handles.SniffTrace.Filtered(pooledsniffs(:,8)));
    set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(pooledsniffs(:,9)),...
        'YData',handles.SniffTrace.Filtered(pooledsniffs(:,9)));

    set(handles.peaksNew, 'XData',[], 'YData',[]);
    set(handles.valleysNew,  'XData',[], 'YData',[]);

    % append timestamps to processed file
    SniffDetectionThreshold = str2double(handles.SDfactor.String);
    CuratedSniffTimestamps = pooledsniffs;

    save(handles.WhereSession.String,'SniffDetectionThreshold','CuratedSniffTimestamps','-append');
    disp('saved sniffs!')
end

function [] = UpdatePeakValleyPlots(hObject, eventdata, handles)
% the original points
olddetections = find(handles.SniffsTS(:,8)>=0);
set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,8)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTS(olddetections,8)));
set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,9)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTS(olddetections,9)));
% new detections
if ~isempty(handles.SniffsTSnew)
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
% hObject    handle to RemovePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi = drawrectangle;
peaks_to_delete = intersect(...
                    find(handles.SniffsTSnew(:,1) >= roi.Position(1)), ...
                        find(handles.SniffsTSnew(:,1) <= (roi.Position(1) + roi.Position(3)) ) );

valleys_to_delete = intersect(...
                    find(handles.SniffsTSnew(:,2) >= roi.Position(1)), ...
                        find(handles.SniffsTSnew(:,2) <= (roi.Position(1) + roi.Position(3)) ) );

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

% --- Executes on button press in ZoomON.
function ZoomON_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom yon
guidata(hObject, handles);

% --- Executes on button press in ZoomOFF.
function ZoomOFF_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomOFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom off
guidata(hObject, handles);

function WindowSize_Callback(hObject, eventdata, handles)
% hObject    handle to WindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowSize as text
%        str2double(get(hObject,'String')) returns contents of WindowSize as a double
currlims = get(handles.SniffingFiltered,'XLim');
newLims = currlims(1) + [0 str2double(handles.WindowSize.String)];
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);
guidata(hObject, handles);


% --- Executes on button press in ClearSession.
function ClearSession_Callback(hObject, eventdata, handles)
% hObject    handle to ClearSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AddSniff.
function AddSniff_Callback(hObject, eventdata, handles)
% hObject    handle to AddSniff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.SniffingFiltered);
zoom off
[x,y] = ginput(2); % select a pair of peak and valley

[~,idx] = min(abs(handles.SniffTrace.Timestamps - x(1)));
% take a window 30 ms on either side
[peakval,peakidx] = max(handles.SniffTrace.Filtered(idx+[-15:1:15]));
handles.newSniff(n,1) = handles.SniffTrace.Timestamps(idx+peakidx-15);

[~,idx] = min(abs(handles.SniffTrace.Timestamps - x(2)));
% take a window 30 ms on either side
[valleyval,valleyidx] = min(handles.SniffTrace.Filtered(idx+[-15:1:15]));
handles.newSniff(n,2) = handles.SniffTrace.Timestamps(idx+valleyidx-15);
guidata(hObject, handles);


% --- Executes on slider movement.
function Scroller_Callback(hObject, eventdata, handles)
% hObject    handle to Scroller (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newLims = handles.Scroller.Value * handles.SessionLength + ...
    [0 str2double(handles.WindowSize.String)];
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in NextStretch.
function NextStretch_Callback(hObject, eventdata, handles)
% hObject    handle to NextStretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newLims = get(handles.SniffingFiltered,'XLim') + 0.9*str2double(handles.WindowSize.String);
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);
handles.Scroller.Value = newLims(1)/(handles.SessionLength - str2double(handles.WindowSize.String));
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in PreviousStretch.
function PreviousStretch_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousStretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newLims = get(handles.SniffingFiltered,'XLim') - 0.9*str2double(handles.WindowSize.String);
set(handles.SniffingRaw,'XLim',newLims);
set(handles.SniffingFiltered,'XLim',newLims);
handles.Scroller.Value = newLims(1)/(handles.SessionLength - str2double(handles.WindowSize.String));
% Update handles structure
guidata(hObject, handles);


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
indices_to_delete = intersect(...
                    find(handles.SniffTrace.Timestamps >= roi.Position(1)), ...
                        find(handles.SniffTrace.Timestamps <= (roi.Position(1) + roi.Position(3)) ) );
handles.rawTrace.YData(indices_to_delete) = nan;
                    
% axes(handles.SniffingRaw);
% set(handles.rawTrace,'XData', handles.SniffTrace.Timestamps, 'YData', handles.SniffTrace.Raw);
delete(roi);
guidata(hObject, handles);


% --- Executes on button press in RedoStretch.
function RedoStretch_Callback(hObject, eventdata, handles)
% hObject    handle to RedoStretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi = drawrectangle;
indices_of_interest = intersect(...
                    find(handles.SniffTrace.Timestamps >= roi.Position(1)), ...
                        find(handles.SniffTrace.Timestamps <= (roi.Position(1) + roi.Position(3)) ) );

sdnew = 2*str2double(handles.SDfactor.String);

switch handles.datamode
    case {'smellocator', 'smellocator_raw'}
        stretchTimeStamps = handles.rawTrace.XData(indices_of_interest)';
        stretchThermistor = handles.rawTrace.YData(indices_of_interest)';
        stretchOdorLocation = handles.OdorLocationTrace(indices_of_interest);

        thisStretchSniffs = ProcessThermistorData([stretchTimeStamps stretchThermistor stretchOdorLocation],'SDfactor',sdnew,'dlgoverride',logical(1));
        thisStretchSniffs = thisStretchSniffs(find(thisStretchSniffs(:,end)>0),:);
        % find trace indices that correspond to detected sniff timestamps
        for n = 1:size(thisStretchSniffs,1)

            % was this already detected
            f = find(abs(handles.SniffsTS(:,1)-thisStretchSniffs(n,1))<0.004,1,'first');
            if ~isempty(f) & abs(handles.SniffsTS(f,2)-thisStretchSniffs(n,2))<0.004
                thisStretchSniffs(n,8:9) = NaN;
            else
                % inhalation start
                [~,idx] = min(abs(handles.SniffTrace.Timestamps - thisStretchSniffs(n,1)));
                if abs(handles.SniffTrace.Timestamps(idx) - thisStretchSniffs(n,1)) < 0.004
                    thisStretchSniffs(n,8) = idx;
                end
                % inhalation end
                [~,idx] = min(abs(handles.SniffTrace.Timestamps - thisStretchSniffs(n,2)));
                if abs(handles.SniffTrace.Timestamps(idx) - thisStretchSniffs(n,2)) < 0.004
                    thisStretchSniffs(n,9) = idx;
                end
            end
        end

        if isfield(handles,'SniffsTSnew')
            handles.SniffsTSnew = vertcat(handles.SniffsTSnew, thisStretchSniffs(find(thisStretchSniffs(:,end)>0),:));
            handles.SniffsTSnew = sortrows(handles.SniffsTSnew,1);
        else
            handles.SniffsTSnew = thisStretchSniffs;
        end

end

% collate the list and make sure peaks fall in order
[handles] = collatesniffs(handles,0);
% update plots
UpdatePeakValleyPlots(hObject, eventdata, handles);

% keep the new detections
answer = questdlg('Keep the new detections?', ...
    'Redo stretch', ...
    'Yes','No','Yes');

% Handle response
switch answer
    case 'Yes'
        
    case 'No'
        for x = 1:size(thisStretchSniffs,1)
            if ~isempty(find(ismember(thisStretchSniffs(x,1:3),handles.SniffsTSnew(:,1:3),'rows')))
                y = find(ismember(thisStretchSniffs(x,1:3),handles.SniffsTSnew(:,1:3),'rows'));
                handles.SniffsTSnew(y,:) = [];
            end
        end
        handles.SniffsTS(find(handles.SniffsTS(:,8)<0),8) = -handles.SniffsTS(find(handles.SniffsTS(:,8)<0),8);
        UpdatePeakValleyPlots(hObject, eventdata, handles); 
end
delete(roi);
guidata(hObject, handles);

