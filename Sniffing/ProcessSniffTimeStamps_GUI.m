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

% Last Modified by GUIDE v2.5 16-Aug-2024 20:49:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcessSniffTimeStamps_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcessSniffTimeStamps_GUI_OutputFcn, ...
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
    if exist(varargin{1}) == 2
        handles.WhereSession.String = varargin{1};
    else
        MouseName = regexprep(varargin{1},'_(\w+)_processed.mat','');
        handles.WhereSession.String = fullfile(Paths.ProcessedSessions,MouseName,varargin{1});
    end
else
    handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
end

% Update handles structure
guidata(hObject, handles);

if exist(handles.WhereSession.String)==2
    LoadSession_Callback(hObject, eventdata, handles);
end

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

% load the relevant data
load(handles.WhereSession.String,'Traces','TrialInfo','SampleRate');
handles.SessionLength = ceil(Traces.Timestamps{end}(end));
handles.SessionDuration.String = num2str(handles.SessionLength);

% reprocess sniff traces on a trialwise basis
Traces.OdorLocation     = Traces.Motor;
[SniffTimeStamps] = ...
    TrialWiseSniffs(TrialInfo,Traces); % [sniffstart sniffstop nextsniff odorlocation sniffslope stimstate trialID]
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

% Trial times
Xvals = [TrialInfo.SessionTimestamps(:,1) TrialInfo.SessionTimestamps(:,1) nan(size(TrialInfo.SessionTimestamps,1),1)]';
YVals = repmat([2*currlims NaN],size(TrialInfo.SessionTimestamps,1),1)';
plot(Xvals(:),YVals(:),'--','color',Plot_Colors('pd'));

Xvals = [TrialInfo.SessionTimestamps(:,2) TrialInfo.SessionTimestamps(:,2) nan(size(TrialInfo.SessionTimestamps,1),1)]';
YVals = repmat([2*currlims NaN],size(TrialInfo.SessionTimestamps,1),1)';
plot(Xvals(:),YVals(:),'-.','color',Plot_Colors('pl'));

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

axes(handles.SniffingFiltered);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in RefindPeaks.
function RefindPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to RefindPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sdnew = 2*str2double(handles.SDfactor.String);

% re-detect peaks and valleys with lower peak prominence
load(handles.WhereSession.String,'Traces','TrialInfo');
Traces.OdorLocation     = Traces.Motor;


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

% % now i can plot all new detections
% newdetections = find(~isnan(handles.SniffsTSnew(:,9)));
% set(handles.peaksNew, ...
%     'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,8)),...
%     'YData',handles.SniffTrace.Filtered(handles.SniffsTSnew(newdetections,8)));
% set(handles.valleysNew, ...
%     'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,9)),...
%     'YData',handles.SniffTrace.Filtered(handles.SniffsTSnew(newdetections,9)));

% collate the list and make sure peaks fall in order
[handles] = collatesniffs(handles,0);
handles = UpdatePeakValleyPlots(handles);
% redraw
% set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),8)),...
%     'YData',handles.SniffTrace.Filtered(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),8)));
% set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),9)),...
%     'YData',handles.SniffTrace.Filtered(handles.SniffsTS(find(handles.SniffsTS(:,8)>=0),9)));

handles.KeepNewSD.BackgroundColor = [0 0.94 0];
handles.KeepOldSD.BackgroundColor = [0 0.94 0];
guidata(hObject, handles);

uiwait;
guidata(hObject, handles);


function [handles] = collatesniffs(handles,whichmode)

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
        f = find(abs(pooledsniffs(1:end-1,3)-pooledsniffs(2:end,1))>0.01);
        for x = 1:numel(f)
            if pooledsniffs(f(x)+1,1)<pooledsniffs(f(x),3) & pooledsniffs(f(x)+1,1)>pooledsniffs(f(x),2) & pooledsniffs(f(x),3)==pooledsniffs(f(x)+1,3)
                pooledsniffs(f(x),3) = pooledsniffs(f(x)+1,1);
            else
                keyboard;
            end
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

function [handles] = UpdatePeakValleyPlots(handles) 
% the original points
olddetections = find(handles.SniffsTS(:,8)>=0);
set(handles.peaksFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,8)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTS(olddetections,8)));
set(handles.valleysFilt,'XData',handles.SniffTrace.Timestamps(handles.SniffsTS(olddetections,9)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTS(olddetections,9)));
% new detections
newdetections = find(handles.SniffsTSnew(:,9)>0);
set(handles.peaksNew, ...
    'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,8)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTSnew(newdetections,8)));
set(handles.valleysNew, ...
    'XData',handles.SniffTrace.Timestamps(handles.SniffsTSnew(newdetections,9)),...
    'YData',handles.SniffTrace.Filtered(handles.SniffsTSnew(newdetections,9)));

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
    handles.SniffsTSnew(peaks_to_delete,9) = -handles.SniffsTSnew(peaks_to_delete,9);
elseif peaks_to_delete == valleys_to_delete + 1
    handles.SniffsTSnew(peaks_to_delete,9) = -handles.SniffsTSnew(peaks_to_delete,9);
    handles.SniffsTSnew(valleys_to_delete,9) = -handles.SniffsTSnew(valleys_to_delete,9);
    % unflag the original sniff in the old detections
    [~,m] = min(abs(handles.SniffsTS(:,1)-handles.SniffsTSnew(valleys_to_delete,1)));
    handles.SniffsTS(m,8) = abs(handles.SniffsTS(m,8));

    handles.lastdeleted.String = mat2str(vertcat(valleys_to_delete,peaks_to_delete,-m));
end
handles = UpdatePeakValleyPlots(handles);    
delete(roi);
guidata(hObject, handles);


% --- Executes on button press in UndoDelete.
function UndoDelete_Callback(hObject, eventdata, handles)
% hObject    handle to UndoDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.lastdeleted.String)
    undelete = eval(handles.lastdeleted.String);
    handles.SniffsTSnew(undelete(undelete>0),9) = abs(handles.SniffsTSnew(undelete(undelete>0),9));
    handles.SniffsTS(abs(undelete(undelete<0)),8) = -abs(handles.SniffsTS(abs(undelete(undelete<0)),8));
    handles = UpdatePeakValleyPlots(handles);
    handles.lastdeleted.String = '';
    guidata(hObject, handles);
end


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
uiresume;
handles.SDfactor.String = num2str(2*str2double(handles.SDfactor.String));
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of KeepNewSD


% --- Executes on button press in KeepOldSD.
function KeepOldSD_Callback(hObject, eventdata, handles)
% hObject    handle to KeepOldSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.peaksNew, 'XData',[], 'YData',[]);
set(handles.valleysNew,  'XData',[], 'YData',[]);
handles.SniffsTS(find(handles.SniffsTS(:,8)<0),8) = -handles.SniffsTS(find(handles.SniffsTS(:,8)<0),8);
handles.KeepNewSD.BackgroundColor = [0.94 0.94 0.94];
handles.KeepOldSD.BackgroundColor = [0.94 0.94 0.94];
uiresume;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of KeepOldSD


% --- Executes on button press in SaveSniffs.
function SaveSniffs_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSniffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles] = collatesniffs(handles,1);
