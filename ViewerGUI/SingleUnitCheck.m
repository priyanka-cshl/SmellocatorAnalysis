function varargout = SingleUnitCheck(varargin)
% SINGLEUNITCHECK MATLAB code for SingleUnitCheck.fig
%      SINGLEUNITCHECK, by itself, creates a new SINGLEUNITCHECK or raises the existing
%      singleton*.
%
%      H = SINGLEUNITCHECK returns the handle to a new SINGLEUNITCHECK or the handle to
%      the existing singleton*.
%
%      SINGLEUNITCHECK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLEUNITCHECK.M with the given input arguments.
%
%      SINGLEUNITCHECK('Property','Value',...) creates a new SINGLEUNITCHECK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SingleUnitCheck_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SingleUnitCheck_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SingleUnitCheck

% Last Modified by GUIDE v2.5 17-Jan-2024 11:51:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SingleUnitCheck_OpeningFcn, ...
                   'gui_OutputFcn',  @SingleUnitCheck_OutputFcn, ...
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


% --- Executes just before SingleUnitCheck is made visible.
function SingleUnitCheck_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SingleUnitCheck (see VARARGIN)

% Choose default command line output for SingleUnitCheck
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

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SingleUnitCheck wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SingleUnitCheck_OutputFcn(hObject, eventdata, handles) 
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
MySession = handles.WhereSession.String;
% Load the relevant variables
load(MySession, 'TrialInfo', 'TTLs', 'ReplayTTLs', 'Tuning*', 'SingleUnits', 'TimestampAdjust', 'SniffTS', 'SniffTS_passive');

% concatenate SniffTS and SniffTS_passive and convert both to OEPS timebase
SniffTS = SniffTS + TimestampAdjust.ClosedLoop;
SniffTS = [SniffTS; (SniffTS_passive + TimestampAdjust.Passive)];

% load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
%     'startoffset', 'errorflags', 'SampleRate', ...
%     'TTLs', 'ReplayTTLs', 'Tuning*', 'SingleUnits');

handles.TrialInfo = TrialInfo;
handles.SingleUnits = SingleUnits;
handles.TimestampAdjuster = TimestampAdjust.ClosedLoop;
handles.SniffTS = SniffTS;

handles.NumUnits.String = num2str(size(SingleUnits,2));
if isnan(handles.CurrentUnit.Data(1)) || handles.CurrentUnit.Data(1)>size(SingleUnits,2)
    handles.CurrentUnit.Data(1) = 1;
end
UpdateUnits(handles);

% Update handles structure
guidata(hObject, handles);

function UpdateUnits(handles)
whichUnit = handles.CurrentUnit.Data(1);

% plot spike amplitudes
axes(handles.SpikeAmplitudes);
thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
hold off
plot(handles.SingleUnits(whichUnit).spikes,thisunitamps,'.','Color',-0.8+[0.8 0.8 0.8],'MarkerSize',10);
hold on
session_end = handles.TrialInfo.SessionTimestamps(end,2) + handles.TimestampAdjuster;
line([session_end session_end],get(gca,'YLim'),'Color','k');

% plot respiration phase
axes(handles.SpikePhases);
thisunitspikes = handles.SingleUnits(whichUnit).spikes;
for i = 1:numel(thisunitspikes)
    thisspike   = thisunitspikes(i);
    whichsniff  = find(handles.SniffTS(:,1)<=thisspike,1,'last') ;
    if ~isempty(whichsniff)
        if thisspike <= handles.SniffTS(whichsniff,3)
            handles.SingleUnits(whichUnit).sniffalignedspikes(i,1) = thisspike -  handles.SniffTS(whichsniff,1); % latency
        else
            handles.SingleUnits(whichUnit).sniffalignedspikes(i,1) = NaN;
        end
    else
        handles.SingleUnits(whichUnit).sniffalignedspikes(i,1) = NaN;
    end
end
thisunitphases = handles.SingleUnits(whichUnit).sniffalignedspikes;
% normalize spike amps 
thisunitamps = thisunitamps - min(thisunitamps,[],'omitnan');
thisunitamps = thisunitamps/max(thisunitamps,[],'omitnan');
hold off
%plot(handles.SingleUnits(whichUnit).spikes,thisunitphases,'.');
colormap('hsv');
scatter(handles.SingleUnits(whichUnit).spikes,thisunitphases,15,thisunitamps,'o','filled');
set(gca,'YLim',eval(handles.RespScaling.String));
hold on
line([session_end session_end],get(gca,'YLim'),'Color','k');


% --- Executes on button press in ClearSession.
function ClearSession_Callback(hObject, eventdata, handles)
% hObject    handle to ClearSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
