function varargout = HaltSniffViewer(varargin)
% HALTSNIFFVIEWER MATLAB code for HaltSniffViewer.fig
%      HALTSNIFFVIEWER, by itself, creates a new HALTSNIFFVIEWER or raises the existing
%      singleton*.
%
%      H = HALTSNIFFVIEWER returns the handle to a new HALTSNIFFVIEWER or the handle to
%      the existing singleton*.
%
%      HALTSNIFFVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HALTSNIFFVIEWER.M with the given input arguments.
%
%      HALTSNIFFVIEWER('Property','Value',...) creates a new HALTSNIFFVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HaltSniffViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HaltSniffViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HaltSniffViewer

% Last Modified by GUIDE v2.5 26-Dec-2023 13:28:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HaltSniffViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @HaltSniffViewer_OutputFcn, ...
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


% --- Executes just before HaltSniffViewer is made visible.
function HaltSniffViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HaltSniffViewer (see VARARGIN)

% Choose default command line output for HaltSniffViewer
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

% UIWAIT makes HaltSniffViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HaltSniffViewer_OutputFcn(hObject, eventdata, handles) 
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

%% check that its a halt session - if not disable Halt only options
if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip')) || ...
            any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
        handles.OnlyHaltRelated.Value = 1; 
        handles.OnlyHaltRelated.Enable = 'on';
        handles.OdorList = mode(...
            [ handles.TrialInfo.Odor(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip')); ...
            handles.TrialInfo.Odor(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))] ...
            );
else
    handles.OnlyHaltRelated.Value = 0; 
    handles.OnlyHaltRelated.Enable = 'off'; 
    handles.OdorList = [1 2 3];
end

%% Get all spikes, all units aligned to trials
[handles.AlignedSniffs, handles.sniffAlignedSpikes, handles.trialAlignedSpikes, ...
    handles.whichtetrode, handles.Events, handles.EventsPhase, handles.TrialInfo] = ...
    TrialAndSniffAlignedSpikeTimes(SingleUnits,TTLs,size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession);

%% same for replays
if ~isempty(OpenLoop)
    
    % sniffs
    [handles.ReplayAlignedSniffs, handles.SniffAlignedReplaySpikes, handles.ReplayInfo] = ...
        SniffAlignedSpikeTimes_Replays(SingleUnits,TTLs,ReplayTTLs,handles.TrialInfo,OpenLoop,MySession);
    
    % regular trials
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
        [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayTrialInfo] = ...
            PerturbationReplayAlignedSpikeTimes_v2(SingleUnits,TTLs,...
            ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop,MySession,'sniffwarpmethod',0);
    end
else
    handles.ReplayAlignedSniffs = [];
end

%% also for passive tuning
handles.TuningSniffs = PassiveTuningSniffs(handles.Tuning,MySession);

%% pseudorandomtuning trials
if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
    [handles.PseudoRandomTuningSpikes] = ...
        TrialAlignedSpikeTimes_Tuning(handles.SingleUnits,handles.Tuning.TTLs);
    
    % transition markers
    odorTS(1,1) = sum(handles.Tuning.extras.sessionsettings(1,4)); % w.r.t. trial start (motor-settle, pre-odor)
    nLocations = size(handles.Tuning.extras.sequence,2) - 2;
    LocationShifts = 0; 
    for i = 1:nLocations
        if i == 1
            LocationShifts(i,1) = -handles.Tuning.extras.sessionsettings(1,3);
            LocationShifts(i,2) = sum(handles.Tuning.extras.sessionsettings(1,[4,5]));
        else
            LocationShifts(i,1) = LocationShifts(i-1,2);
            LocationShifts(i,2) = LocationShifts(i,1) + sum(handles.Tuning.extras.sessionsettings(1,[4,5]));
        end
    end
    odorTS(1,2) = LocationShifts(end,2);
    handles.TuningTiming.LocationShifts = LocationShifts/1000; % in s
    handles.TuningTiming.Odor = odorTS/1000; % in s

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
MyColors1 = brewermap(15,'*PuBu');
MyColors2 = brewermap(15,'*OrRd');
handles.tetrode.String = num2str(handles.whichtetrode(whichUnit,1));
handles.Cluster_ID.String = num2str(handles.whichtetrode(whichUnit,2));

myXlim = eval(handles.xlims.String); %[-0.1 1.1];

for i = 1:numel(handles.OdorList)
    axes(handles.(['axes',num2str(i)])); 
    cla reset; 
    hold on
    
    whichodor = handles.OdorList(i);
    
    % plot baseline trials
    [nSniffs,~,~,haltlocation] = PlotSortedSniffs(whichUnit, whichodor, handles.trialAlignedSpikes, handles.AlignedSniffs, ...
                                 handles.TrialInfo, 'plotspikes', 0, 'sortorder', (handles.SortOrder.Value-1), ...
                                 'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1, ...
                                 'haltmode', handles.OnlyHaltRelated.Value);
                             
    
    % plot passive replay trials                         
    if ~isempty(handles.ReplayAlignedSniffs)
        [nSniffs] = PlotPassiveReplaySniffs(nSniffs, whichUnit, whichodor, handles.SniffAlignedReplaySpikes, handles.ReplayAlignedSniffs, ...
            handles.ReplayInfo, 'plotspikes', 0, 'sortorder', (handles.SortOrder.Value-1), ...
            'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1, ...
                                 'haltmode', handles.OnlyHaltRelated.Value);
    end
    
    % add tuning sniffs
    [nSniffs] = PlotTuningSniffs(whichUnit, whichodor, handles.SingleUnits, handles.TuningSniffs, handles.Tuning.extras.sequence, nSniffs, ...
        'plotspikes', 0, 'sortorder', (handles.SortOrder.Value-1), ...
        'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1, ...
        'selectlocation', haltlocation);
            
    if nSniffs>0
        set(gca, 'XLim', myXlim, 'YLim', [0 nSniffs], 'YTick', []);
    end
end

if handles.OnlyHaltRelated.Value
    % add the trial aligned plot
    AlignType = 6; % perturbation start
    myXlim = [-1.2 6];
    
    axes(handles.(['axes',num2str(i+1)])); 
    cla reset; 
    hold on
    % plot baseline trials
    [trialsdone] = PlotFullSession(whichUnit, whichodor, handles.trialAlignedSpikes, handles.Events, ...
        handles.TrialInfo, handles.TrialInfo.InZone, AlignType, 'plotspikes', 0, ...
        'trialfilter', handles.PlotSelectTrials.Value);
    
    % passive halts
    if ~isempty(handles.ReplayAlignedSniffs)
        [perturbationreplaysadded] = AddPerturbationReplay2FullSession_v2(trialsdone, whichUnit, whichodor, handles.ReplayAlignedSpikes, ...
            handles.ReplayEvents, handles.ReplayTrialInfo, handles.ReplayTrialInfo.InZone, AlignType, handles.SortReplay.Value, ...
            'trialfilter', handles.PlotSelectTrials.Value, 'plotspikes', 0);
        
        trialsdone = trialsdone + perturbationreplaysadded;
    end
    
    if handles.PlotCenterTuning.Value
        haltlocation = 0;
    end
    % add tuning trials
    if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
        AlignType = 1000 + haltlocation;
        LocationDuration = mode(diff(handles.TuningTiming.LocationShifts'));
        [trialsdone] = PlotRandomTuningTrials(trialsdone, whichUnit, whichodor, handles.PseudoRandomTuningSpikes, ...
            handles.TuningTiming, handles.Tuning.extras.sequence, AlignType, LocationDuration, myXlim, 'plotspikes', 0);
    else
        [trialsdone] = PlotTuningTrials(trialsdone, whichUnit, whichodor, handles.SingleUnits, handles.Tuning.TTLs, ...
            'plotspikes', 0, 'selectlocation', haltlocation);
    end
 
    set(gca, 'XLim', myXlim, 'YLim', [0 (trialsdone)]);
    %set(gca, 'XLim', myXlim, 'YLim', [0 (trialsdone + replaysadded + perturbationreplaysadded)]);
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
for i = 1:numel(handles.OdorList)
    axes(handles.(['axes',num2str(i+9)])); 
    cla reset; 
    set(gca,'color','none');
    hold on
    
    whichodor = handles.OdorList(i);
    
    % plot baseline trials
    [nSniffs,FR,BinOffset,haltlocation] = PlotSortedSniffs(whichUnit, whichodor, handles.trialAlignedSpikes, handles.AlignedSniffs, ...
                     handles.TrialInfo, 'plotevents', 0, 'sortorder', (handles.SortOrder.Value-1), ...
                     'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1, ...
                     'psth', handles.plotPSTH.Value, 'haltmode', handles.OnlyHaltRelated.Value);
                 
    
    % plot passive replay trials                         
    if ~isempty(handles.ReplayAlignedSniffs)
        [nSniffs,PR_FR] = PlotPassiveReplaySniffs(nSniffs, whichUnit, whichodor, handles.SniffAlignedReplaySpikes, handles.ReplayAlignedSniffs, ...
            handles.ReplayInfo, 'plotevents', 0, 'sortorder', (handles.SortOrder.Value-1), ...
            'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1, ...
            'psth', handles.plotPSTH.Value, 'haltmode', handles.OnlyHaltRelated.Value);
    end
    
    % add tuning sniffs
    [~,T_FR] = PlotTuningSniffs(whichUnit, whichodor, handles.SingleUnits, handles.TuningSniffs, handles.Tuning.extras.sequence, nSniffs, ...
        'plotevents', 0, 'sortorder', (handles.SortOrder.Value-1), ...
        'alignto', handles.SniffAlignment.Value, 'warptype', handles.WarpType.Value-1, ...
        'psth', handles.plotPSTH.Value, 'selectlocation', haltlocation);
    
    set(gca, 'XLim', myXlim, 'YTick', []);
    set(gca, 'YLim', handles.(['axes',num2str(i)]).YLim);
   
    % plot PSTHs
    if handles.plotPSTH.Value
        axes(handles.(['axes',num2str(i+3)]));
        cla reset;
        hold on
        
        if handles.OnlyHaltRelated.Value
            % plot the close loop sniffs
            plot((1:size(FR{5},1))*0.002+BinOffset/1000,FR{5},'Linewidth',2,'Color','k');
            % plot active halt sniffs
            plot((1:size(FR{6},1))*0.002+BinOffset/1000,FR{6},'Linewidth',2,'Color',Plot_Colors('t'));
            if ~isempty(handles.ReplayAlignedSniffs)
                % plot the passive replay control sniffs
                plot((1:size(PR_FR{5},1))*0.002+BinOffset/1000,PR_FR{5},'Linewidth',2,'Color',Plot_Colors('b'));
                % plot the passive halt sniffs
                plot((1:size(PR_FR{6},1))*0.002+BinOffset/1000,PR_FR{6},'Linewidth',2,'Color',Plot_Colors('r'));
            end
            % plot tuning sniffs
            plot((1:size(T_FR{whichodor+2},1))*0.002+BinOffset/1000,T_FR{whichodor+2},'Linewidth',2,'Color',Plot_Colors('o'));
        else
            
            for t = 1:size(FR,2)
                %set(groot,'defaultAxesColorOrder',MyColors2);
                if t < 5
                    plot((1:size(FR{t},1))*0.002+BinOffset/1000,FR{t},'Linewidth',2,'Color',MyColors1(t+3,:));
                else
                    plot((1:size(FR{t},1))*0.002+BinOffset/1000,FR{t},'Linewidth',2,'Color','k');
                end
            end
        end
        set(gca, 'XLim', myXlim);
        myYlims(i,:) = get(gca, 'YLim');
    end
end

% if handles.plotPSTH.Value
%     for i = 1:3
%         set(handles.(['axes',num2str(i+3)]),'YLim', [min(myYlims(:,1)) max(myYlims(:,2))]);
%     end
% end

if handles.OnlyHaltRelated.Value
    % add the trial aligned plot
    AlignType = 6; % perturbation start
    %myXlim = [-1.2 2];
    myXlim = [-1.2 6];
    
    axes(handles.(['axes',num2str(i+9+1)])); 
    cla reset; 
    set(gca,'color','none');
    hold on
    
    % plot baseline trials
    [trialsdone, FRs, BinOffset, P_FRs] = PlotFullSession(whichUnit, whichodor, handles.trialAlignedSpikes, handles.Events, ...
        handles.TrialInfo, handles.TrialInfo.InZone, AlignType, 'plotevents', 0, ...
        'trialfilter', handles.PlotSelectTrials.Value, 'psth', handles.plotPSTH.Value, 'poolTZs', handles.poolTZs.Value);
    
    % passive halts
    if ~isempty(handles.ReplayAlignedSniffs)
        [perturbationreplaysadded, replayFRs, perturbreplayFRs, replayOffset] = AddPerturbationReplay2FullSession_v2(trialsdone, whichUnit, whichodor, handles.ReplayAlignedSpikes, ...
            handles.ReplayEvents, handles.ReplayTrialInfo, handles.ReplayTrialInfo.InZone, AlignType, handles.SortReplay.Value, ...
            'trialfilter', handles.PlotSelectTrials.Value, 'plotevents', 0, 'psth', handles.plotPSTH.Value, 'poolTZs', handles.poolTZs.Value);
    else
        perturbationreplaysadded = 0;
    end
    trialsdone = trialsdone + perturbationreplaysadded;
    
    if handles.PlotCenterTuning.Value
        haltlocation = 0;
    end
    % add tuning trials
    if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
        AlignType = 1000 + haltlocation;
        LocationDuration = mode(diff(handles.TuningTiming.LocationShifts'));
        [trialsdone, TuningFR, TuningOffset] = PlotRandomTuningTrials(trialsdone, whichUnit, whichodor, handles.PseudoRandomTuningSpikes, ...
            handles.TuningTiming, handles.Tuning.extras.sequence, AlignType, LocationDuration, myXlim, 'plotevents', 0, 'psth', handles.plotPSTH.Value);
    else
        [trialsdone, TuningFR, TuningOffset] = PlotTuningTrials(trialsdone, whichUnit, whichodor, handles.SingleUnits, handles.Tuning.TTLs, ...
        'plotevents', 0, 'selectlocation', haltlocation, 'psth', handles.plotPSTH.Value);
    end
    %TuningFR = [];
    
    set(gca, 'XLim', myXlim);
    set(gca, 'YLim', handles.(['axes',num2str(i+1)]).YLim);
end

if handles.plotPSTH.Value

    axes(handles.(['axes',num2str(i+3+1)]));
    cla reset;
    hold on
    
    if ~any(strcmp(handles.TrialInfo.Perturbation(:,1),'RuleReversal'))
        for t = 1:size(FRs,1)
            plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(t,:),'Color','k','Linewidth',2);
        end
        
        if ~isempty(P_FRs)
            for t = 1:size(P_FRs,1)
                set(groot,'defaultAxesColorOrder',MyColors2);
                plot((1:size(P_FRs,2))*0.002+BinOffset/1000,P_FRs(t,:),'Color',Plot_Colors('t'),'Linewidth',2);
            end
        end
        
        % passivereplays and passive halts
        if perturbationreplaysadded
            if ~isempty(replayFRs)
                for t = 1:size(replayFRs,2)
                    plot((1:numel(replayFRs{t}))*0.002+replayOffset/1000,replayFRs{t},'Color',Plot_Colors('b'),'Linewidth',2);
                end
            end
            
            if ~isempty(perturbreplayFRs)
                for t = 1:size(perturbreplayFRs,2)
                    plot((1:numel(perturbreplayFRs{t}))*0.002+replayOffset/1000,perturbreplayFRs{t},'Color',Plot_Colors('r'),'Linewidth',2);
                end
            end
        end
    
        if ~isempty(TuningFR)
            plot((1:size(TuningFR,2))*0.002+TuningOffset/1000,TuningFR(1,:),'Color',Plot_Colors('o'),'Linewidth',2);
        end
        
    end
    set(gca, 'XLim', myXlim);
end
    
% plot spike amplitudes
if handles.spike_amplitudes.Value
    thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
    axes(handles.amplitudeaxes);
    hold off
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


% --- Executes on button press in SortReplay.
function SortReplay_Callback(hObject, eventdata, handles)
% hObject    handle to SortReplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SortReplay
