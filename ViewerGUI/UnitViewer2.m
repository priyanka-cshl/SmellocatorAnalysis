function varargout = UnitViewer2(varargin)
% UNITVIEWER2 MATLAB code for UnitViewer2.fig
%      UNITVIEWER2, by itself, creates a new UNITVIEWER2 or raises the existing
%      singleton*.
%
%      H = UNITVIEWER2 returns the handle to a new UNITVIEWER2 or the handle to
%      the existing singleton*.
%
%      UNITVIEWER2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNITVIEWER2.M with the given input arguments.
%
%      UNITVIEWER2('Property','Value',...) creates a new UNITVIEWER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UnitViewer2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UnitViewer2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UnitViewer2

% Last Modified by GUIDE v2.5 18-Aug-2023 11:27:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UnitViewer2_OpeningFcn, ...
                   'gui_OutputFcn',  @UnitViewer2_OutputFcn, ...
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


% --- Executes just before UnitViewer2 is made visible.
function UnitViewer2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UnitViewer2 (see VARARGIN)

% Choose default command line output for UnitViewer2
handles.output = hObject;

[Paths] = WhichComputer();

if ~isempty(varargin)
    MouseName = regexprep(varargin{1},'_(\w+)_processed.mat','');
    handles.WhereSession.String = fullfile(Paths.ProcessedSessions,MouseName,varargin{1});
else
    handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
end

handles.CurrentUnit.Data(1) = NaN;

set(handles.axes10,'Color','none');
set(handles.axes11,'Color','none');
set(handles.axes12,'Color','none');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UnitViewer2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UnitViewer2_OutputFcn(hObject, eventdata, handles) 
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
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

handles.SingleUnits = SingleUnits;
handles.TimestampAdjuster = TimestampAdjuster;

%% Get all spikes, all units aligned to trials
[handles.AlignedSpikes, handles.Events, handles.whichtetrode, handles.sniffAlignedSpikes, handles.EventsPhase, handles.TrialInfo] = ...
    TrialAlignedSpikeTimes_v2(SingleUnits,TTLs,...
    size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession,'sniffwarpmethod',handles.SniffAlignment.Value-1);

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo, handles.SniffAlignedReplaySpikes, handles.SniffAlignedReplayEvents] = ...
        PerturbationReplayAlignedSpikeTimes_v2(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop,MySession,'sniffwarpmethod',handles.SniffAlignment.Value-1);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Offset-II-Template'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop);
end

handles.NumUnits.String = num2str(size(SingleUnits,2));
if isnan(handles.CurrentUnit.Data(1)) || handles.CurrentUnit.Data(1)>size(SingleUnits,2)
    handles.CurrentUnit.Data(1) = 1;
end

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

if handles.AlignToSniffs.Value
    myXlim = str2num(handles.sniffscalar.String)*myXlim;
end

for i = 1:3
    axes(handles.(['axes',num2str(i)])); 
    cla reset; 
    hold on
    % plot baseline trials
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'RuleReversal'))
        [trialsdone, FRs, BinOffset, P_FRs] = PlotRuleReversalSession(whichUnit, i, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignType);
    else
        if handles.AlignToSniffs.Value
            [trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, i, handles.sniffAlignedSpikes, handles.EventsPhase, ...
                handles.TrialInfo, handles.TrialInfo.InZonePhase, AlignType, 'plotspikes', 0, 'sniffaligned', 1, 'sniffscalar', str2num(handles.sniffscalar.String));
        else
            [trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, i, handles.AlignedSpikes, handles.Events, ...
                handles.TrialInfo, handles.TrialInfo.InZone, AlignType, 'plotspikes', 0);
        end
    end
    if any(strcmp(handles.TrialInfo.Perturbation,'OL-Replay'))
        % plot replay trials
        AddReplay2FullSession(trialsdone, whichUnit, i, handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo, AlignType, handles.SortReplay.Value, 'plotspikes', 0);
    end
    
    % plot halt/offset replays
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template')) || ...
            any(strcmp(handles.TrialInfo.Perturbation(:,1),'Offset-II-Template'))
        % plot passive halt trials
        if handles.AlignToSniffs.Value
            AddPerturbationReplay2FullSession(trialsdone, whichUnit, i, handles.SniffAlignedReplaySpikes, ...
                handles.SniffAlignedReplayEvents, handles.ReplayInfo, handles.ReplayInfo.InZonePhase, AlignType, handles.SortReplay.Value,...
                'plotspikes', 0, 'sniffaligned', 1, 'sniffscalar', str2num(handles.sniffscalar.String));
        else
            AddPerturbationReplay2FullSession(trialsdone, whichUnit, i, handles.ReplayAlignedSpikes, ...
                handles.ReplayEvents, handles.ReplayInfo, handles.ReplayInfo.InZone, AlignType, handles.SortReplay.Value, 'plotspikes', 0);
        end
    end
    
    set(gca, 'XLim', myXlim);
    
%     axes(handles.(['axes',num2str(i+3)])); 
%     cla reset; 
%     hold on
%     
%     if ~any(strcmp(handles.TrialInfo.Perturbation(:,1),'RuleReversal'))
%         for t = 1:size(FRs,1)
%             plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(t,:),'Color',MyColors1(t,:),'Linewidth',1);
%         end
%         
%         if ~handles.HidePSTH2.Value
%             if ~isempty(P_FRs)
%                 for t = 1:size(FRs,1)
%                     set(groot,'defaultAxesColorOrder',MyColors2);
%                     plot((1:size(P_FRs,2))*0.002+BinOffset/1000,P_FRs(t,:),'Color',MyColors2(t,:),'Linewidth',1);
%                 end
%             end
%         end
%     end
%     set(gca, 'XLim', myXlim);
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
if handles.AlignToSniffs.Value
    myXlim = str2num(handles.sniffscalar.String)*myXlim;
end
for i = 1:3
    axes(handles.(['axes',num2str(i+9)])); 
    cla reset; 
    set(gca,'color','none');
    hold on
    % plot baseline trials
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'RuleReversal'))
        [trialsdone, FRs, BinOffset, P_FRs] = PlotRuleReversalSession(whichUnit, i, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignType);
    else
        if handles.AlignToSniffs.Value
            [trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, i, handles.sniffAlignedSpikes, handles.EventsPhase, ...
                handles.TrialInfo, handles.TrialInfo.InZonePhase, AlignType, 'plotevents', 0, 'sniffaligned', 1, 'sniffscalar', str2num(handles.sniffscalar.String));
        else
            [trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, i, handles.AlignedSpikes, handles.Events, ...
                handles.TrialInfo, handles.TrialInfo.InZone, AlignType, 'plotevents', 0);
        end
    end
    if any(strcmp(handles.TrialInfo.Perturbation,'OL-Replay'))
        % plot replay trials
        AddReplay2FullSession(trialsdone, whichUnit, i, handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo, AlignType, handles.SortReplay.Value, 'plotevents', 0);
    end
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template')) || ...
            any(strcmp(handles.TrialInfo.Perturbation(:,1),'Offset-II-Template'))
         % plot passive halt trials
        if handles.AlignToSniffs.Value
            AddPerturbationReplay2FullSession(trialsdone, whichUnit, i, handles.SniffAlignedReplaySpikes, ...
                handles.SniffAlignedReplayEvents, handles.ReplayInfo, handles.ReplayInfo.InZonePhase, AlignType, handles.SortReplay.Value,...
                'plotevents', 0, 'sniffaligned', 1, 'sniffscalar', str2num(handles.sniffscalar.String));
        else
            AddPerturbationReplay2FullSession(trialsdone, whichUnit, i, handles.ReplayAlignedSpikes, ...
                handles.ReplayEvents, handles.ReplayInfo, handles.ReplayInfo.InZone, AlignType, handles.SortReplay.Value, 'plotevents', 0);
        end
    end
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
    
    if handles.AlignToSniffs.Value
        % sniff phase related analysis
        axes(handles.(['axes',num2str(i+13)]));
        hold off
        PlotPreferredPhase(whichUnit, i, handles.sniffAlignedSpikes, handles.EventsPhase, handles.TrialInfo);
    end
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


% --- Executes on button press in plotPSTH.
function plotPSTH_Callback(hObject, eventdata, handles)
% hObject    handle to plotPSTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotPSTH
