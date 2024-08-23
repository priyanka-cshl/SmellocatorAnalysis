function varargout = Smellocator_DataViewer_PG(varargin)
% SMELLOCATOR_DATAVIEWER_PG MATLAB code for Smellocator_DataViewer_PG.fig
%      SMELLOCATOR_DATAVIEWER_PG, by itself, creates a new SMELLOCATOR_DATAVIEWER_PG or raises the existing
%      singleton*.
%
%      H = SMELLOCATOR_DATAVIEWER_PG returns the handle to a new SMELLOCATOR_DATAVIEWER_PG or the handle to
%      the existing singleton*.
%
%      SMELLOCATOR_DATAVIEWER_PG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMELLOCATOR_DATAVIEWER_PG.M with the given input arguments.
%
%      SMELLOCATOR_DATAVIEWER_PG('Property','Value',...) creates a new SMELLOCATOR_DATAVIEWER_PG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Smellocator_DataViewer_PG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Smellocator_DataViewer_PG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Smellocator_DataViewer_PG

% Last Modified by GUIDE v2.5 10-Aug-2023 09:58:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Smellocator_DataViewer_PG_OpeningFcn, ...
    'gui_OutputFcn',  @Smellocator_DataViewer_PG_OutputFcn, ...
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


% --- Executes just before Smellocator_DataViewer_PG is made visible.
function Smellocator_DataViewer_PG_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% defaults
handles.SampleRate = 500;
handles.SessionLength.String = '100';
handles.TimeWindow.String = '20';
handles.RespirationScaling.Data = [6 0.5];
if ~isempty(varargin)
    handles.WhereSession.String = varargin{1};
else
    handles.WhereSession.String = '';
end
handles.Scroller.Value = 0.2;

% define plot colors
handles.plotcolors.Odor1    = [.8 .8 .8];
handles.plotcolors.Odor2    = [0.8941    0.9412    0.9020];
handles.plotcolors.Odor3    = [0.8706    0.9216    0.9804];
handles.plotcolors.Odor4    = [0.93    0.84    0.84];
handles.plotcolors.TZ       = [1 1 0];
handles.plotcolors.resp     = [52 101 164]./256;
handles.plotcolors.licks    = [239 41 41]./256;
handles.plotcolors.rewards  = [0 139 139]./256;
handles.OdorLocationMapping = [ ...
    0.6471    0.0000    0.1490
    0.6707    0.0081    0.1416
    0.6934    0.0211    0.1361
    0.7152    0.0398    0.1324
    0.7362    0.0612    0.1304
    0.7563    0.0823    0.1303
    0.7756    0.1035    0.1319
    0.7941    0.1249    0.1351
    0.8118    0.1465    0.1399
    0.8287    0.1683    0.1462
    0.8447    0.1905    0.1538
    0.8599    0.2130    0.1625
    0.8743    0.2358    0.1722
    0.8878    0.2590    0.1827
    0.9004    0.2826    0.1939
    0.9123    0.3066    0.2055
    0.9232    0.3310    0.2175
    0.9333    0.3558    0.2296
    0.9425    0.3810    0.2417
    0.9508    0.4066    0.2535
    0.9583    0.4327    0.2650
    0.9649    0.4592    0.2760
    0.9706    0.4858    0.2868
    0.9755    0.5125    0.2974
    0.9797    0.5391    0.3081
    0.9832    0.5655    0.3191
    0.9861    0.5914    0.3305
    0.9883    0.6169    0.3427
    0.9901    0.6418    0.3557
    0.9914    0.6659    0.3698
    0.9924    0.6893    0.3851
    0.9931    0.7117    0.4018
    0.9936    0.7333    0.4195
    0.9938    0.7543    0.4381
    0.9940    0.7745    0.4572
    0.9941    0.7942    0.4767
    0.9942    0.8133    0.4964
    0.9944    0.8319    0.5159
    0.9948    0.8501    0.5352
    0.9955    0.8679    0.5538
    0.9966    0.8854    0.5718
    0.9980    0.9025    0.5890
    0.9996    0.9189    0.6057
    1.0000    0.9345    0.6223
    1.0000    0.9491    0.6391
    1.0000    0.9623    0.6565
    1.0000    0.9741    0.6746
    1.0000    0.9841    0.6939
    1.0000    0.9921    0.7146
    1.0000    0.9980    0.7371
    0.9979    1.0000    0.7615
    0.9924    1.0000    0.7876
    0.9851    1.0000    0.8146
    0.9760    0.9984    0.8418
    0.9651    0.9938    0.8685
    0.9524    0.9880    0.8939
    0.9379    0.9810    0.9172
    0.9220    0.9733    0.9379
    0.9046    0.9650    0.9553
    0.8861    0.9564    0.9685
    0.8667    0.9477    0.9771
    0.8468    0.9390    0.9810
    0.8263    0.9302    0.9810
    0.8055    0.9211    0.9774
    0.7845    0.9116    0.9710
    0.7633    0.9016    0.9624
    0.7419    0.8911    0.9521
    0.7204    0.8798    0.9408
    0.6988    0.8678    0.9289
    0.6771    0.8550    0.9172
    0.6553    0.8412    0.9060
    0.6335    0.8266    0.8955
    0.6116    0.8110    0.8856
    0.5897    0.7947    0.8762
    0.5678    0.7776    0.8671
    0.5459    0.7597    0.8581
    0.5241    0.7413    0.8492
    0.5023    0.7222    0.8403
    0.4806    0.7026    0.8311
    0.4592    0.6825    0.8216
    0.4379    0.6619    0.8116
    0.4170    0.6410    0.8012
    0.3965    0.6196    0.7904
    0.3765    0.5979    0.7793
    0.3570    0.5758    0.7678
    0.3383    0.5534    0.7561
    0.3203    0.5307    0.7441
    0.3032    0.5078    0.7319
    0.2871    0.4846    0.7196
    0.2720    0.4612    0.7071
    0.2581    0.4375    0.6946
    0.2455    0.4137    0.6820
    0.2341    0.3896    0.6695
    0.2241    0.3653    0.6569
    0.2155    0.3408    0.6444
    0.2082    0.3160    0.6321
    0.2023    0.2908    0.6198
    0.1977    0.2651    0.6078
    0.1944    0.2389    0.5959
    0.1922    0.2118    0.5843 ];

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Smellocator_DataViewer_PG_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)

global TargetZones;
global SampleRate;
global startoffset;

if isempty(handles.WhereSession.String)
    disp('Please specify the behavior session to plot');
    return;
end

% Load the relevant variables
load(handles.WhereSession.String, 'Traces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'SampleRate');

handles.SampleRate = SampleRate;

%%
% create a Trace which simply gates periods by sniff Timestamps
if ~isempty(handles.RespSegmentsPath.String)
    sniffTS = load(handles.RespSegmentsPath.String);
    lag1 = round(str2double(handles.lag1.String)*SampleRate);
    lag2 = round(str2double(handles.lag2.String)*SampleRate);
    
    for trial = 1:size(TrialInfo.TrialID,1)
        ts = [];
        x = sniffTS.(handles.phase1.String).(['n',num2str(trial-1)]);
        y = sniffTS.(handles.phase2.String).(['n',num2str(trial-1)]);
        for phases = 1:numel(x)-1
            f = intersect(find(y>x(phases)),find(y<x(phases+1)));
            if ~isempty(f)
                ts(phases,2) = y(f(1));
                if numel(f)>1
                    disp(trial-1);
                end
            else
                ts(phases,2) = x(phases);
            end
            ts(phases,1) = x(phases);
        end
        ts(find((ts(:,1)-ts(:,2))==0),:) = [];
        ts(:,3) = ts(:,1) + lag1;
        ts(:,4) = ts(:,2) + lag2;
        Traces.SniffGated{trial} = 0*Traces.Lever{trial};
        Traces.SniffDelayed{trial} = 0*Traces.Lever{trial};
        for phases = 1:size(ts,1)
            Traces.SniffGated{trial}(ts(phases,1):ts(phases,2)) = 1;
            Traces.SniffDelayed{trial}(ts(phases,3):ts(phases,4)) = 1;
        end
    end
    
end

%%
[TracesOut, whichTraces] = ConcatenateTraces2Matrix(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);
% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(whichTraces,'Timestamps')));

% Sniffing specifc
% add a filtered sniff trace
TracesOut(:,9)     = FilterThermistor(TracesOut(:,2));
whichTraces{9,1}   = 'SniffsFiltered';

% add a digital sniff trace
load(handles.WhereSession.String,'CuratedSniffTimestamps');
if exist('CuratedSniffTimestamps','var')
    if size(CuratedSniffTimestamps,2) < 10
        CuratedSniffTimestamps(:,10) = 0;
    end
    LocationSniffs = TracesOut(:,9)*nan;
    DigitalSniffs = TracesOut(:,9)*0;
    for n = 1:size(CuratedSniffTimestamps)
        idx = CuratedSniffTimestamps(n,8:9);
        if CuratedSniffTimestamps(n,10) == -1
            DigitalSniffs(idx(1):idx(2)) = -1;
        else
            DigitalSniffs(idx(1):idx(2)) = 1;
        end
        if ~any(isnan(CuratedSniffTimestamps(:,4)))
            location = CuratedSniffTimestamps(n,4);
            LocationSniffs(idx(1):idx(2)) = location;
        end
    end

    TracesOut(:,10)     = DigitalSniffs;
    whichTraces{10,1}   = 'SniffsDigitized';
    TracesOut(:,11)     = LocationSniffs;
    whichTraces{11,1}   = 'SniffsLocationed';
    TracesOut(:,12)     = TracesOut(:,1); % Lever
    TracesOut(find(~TracesOut(:,10)),12) = nan;
end

%% Calculate Trial On-Off timestamps
TrialColumn = TracesOut(:,find(strcmp(whichTraces,'TrialState')));
TrialColumn(TrialColumn~=0) = 1; % make logical
TrialOn = find(diff([0; TrialColumn])>0);
TrialOff =  find(diff(TrialColumn)<0)+1;

% account for cases where acquisition started/ended in between a trial
while TrialOn(1)>TrialOff(1)
    TrialOff(1,:) = [];
end
while TrialOn(end)>TrialOff(end)
    TrialOn(end,:) = [];
end

TrialIndices = [TrialOn TrialOff];
TrialTimeStamps = [Timestamps(TrialOn,1) Timestamps(TrialOff,1)];

handles.SessionLength.String = TrialTimeStamps(end,2);

%% plot odor boxes on the behavior plot
axes(handles.BehaviorPlot);
hold off
for i = 1:4
    handles.(['Trial',num2str(i),'Plot']) = fill(NaN,NaN,handles.plotcolors.(['Odor',num2str(i)]));
    hold on;
    handles.(['Trial',num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TrialTimeStamps((TrialInfo.Odor==i),1:2)';
    if ~isempty(ValveTS)
        handles.(['Trial',num2str(i),'Plot']).Vertices = [ ...
            reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
            repmat([0 5 5 0]',size(ValveTS,2),1)];
        handles.(['Trial',num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
    end
end

% plot the target zone
handles.TargetZonePlot = fill(NaN,NaN,handles.plotcolors.TZ,'FaceAlpha',0.2);
hold on;
handles.TargetZonePlot.EdgeColor = 'none';
TrialTS = TrialTimeStamps(:,1:2)';
TZList =  TargetZones(TrialInfo.TargetZoneType,[3 1 1 3])';
handles.TargetZonePlot.Vertices = [ ...
    reshape([TrialTS(:) TrialTS(:)]', 2*numel(TrialTS), []) , ...
    TZList(:)];
handles.TargetZonePlot.Faces = reshape(1:2*numel(TrialTS),4,size(TrialTS,2))';

% plot the lever trace on top
plot(Timestamps, TracesOut(:,find(strcmp(whichTraces,'Lever'))),'k','Linewidth',2);

% plot Rewards
handles.reward_plot = plot(NaN, NaN, 'color',handles.plotcolors.rewards,'Linewidth',1.25);
tick_timestamps =  Timestamps(find(diff(TracesOut(:,find(strcmp(whichTraces,'Rewards')))==1)) + 1);
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; 5; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.reward_plot,'XData',tick_x,'YData',tick_y);

% plot Licks
handles.lick_plot = plot(NaN, NaN, 'color',handles.plotcolors.licks,'Linewidth',0.25);
if ~isempty(find(strcmp(whichTraces,'LicksBinary')))
    tick_timestamps =  Timestamps(find(diff(TracesOut(:,find(strcmp(whichTraces,'LicksBinary')))==1)) + 1);
elseif ~isempty(find(strcmp(whichTraces,'Licks')))
    tick_timestamps =  Timestamps(find(diff(TracesOut(:,find(strcmp(whichTraces,'Licks')))==1)) + 1);
end
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; 5.5; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.lick_plot,'XData',tick_x,'YData',tick_y);

set(gca,'YLim', [0 7], 'YTick', [],...
    'XTick', [], 'XLim', [0 str2double(handles.TimeWindow.String)]);

% plot respiration
handles.SniffData = TracesOut(:,find(strcmp(whichTraces,'Sniffs')));
handles.SniffData(handles.SniffData==Inf) = NaN;
handles.SniffData = handles.SniffData - mean(handles.SniffData,'omitnan');

handles.respiration_plot = plot(Timestamps, ...
    handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.SniffData, ...
    'color',handles.plotcolors.resp,'Linewidth',2);

%%

if ~isempty(handles.RespSegmentsPath.String)
    % gated respiration plot
    handles.gatedsniffs = handles.SniffData;
    handles.gatedsniffs(find(TracesOut(:,9)~=1)) = NaN;
    handles.respiration_plot_gated = plot(Timestamps, ...
        handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.gatedsniffs, ...
        'color',handles.plotcolors.licks,'Linewidth',2);
    
    % gated Lever plot
    gatedLever = TracesOut(:,find(strcmp(whichTraces,'Lever')));
    gatedLever(find(TracesOut(:,10)~=1)) = NaN;
    handles.lever_plot_gated = plot(Timestamps, ...
        gatedLever, ...
        'color',handles.plotcolors.licks,'Linewidth',2);
end

if size(TracesOut,2) == 12
        handles.lever_plot_gated = plot(Timestamps, ...
        TracesOut(:,12), ...
        'color',handles.plotcolors.licks,'Linewidth',2);
end

%% update the motor plot
axes(handles.MotorPlot);
hold off
MotorTrajectory = NaN + zeros(str2double(handles.SessionLength.String)*SampleRate,1);
% Fill in the Motor Data from the behavior session period
Indices = round((Timestamps)*SampleRate);
if ~Indices(1)
    Indices = Indices + 1;
end
MotorTrajectory(Indices,1) = TracesOut(:,find(strcmp(whichTraces,'OdorLocation')))/100;

alphamask = ~isnan(MotorTrajectory);

handles.MotorTrajectoryPlot = imagesc(MotorTrajectory',[-1 1]);
colormap(handles.OdorLocationMapping);
set(handles.MotorTrajectoryPlot, 'AlphaData', alphamask');

set(gca, 'YTick', [], 'XTick', [], ...
    'TickDir','out','XLim', [0 SampleRate*str2double(handles.TimeWindow.String)]);

Scroller_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in NewSession.
function NewSession_Callback(hObject, eventdata, handles)
handles.WhereSession.String = '';

% --- Executes on slider movement.
function Scroller_Callback(hObject, eventdata, handles)
newLims = handles.Scroller.Value * str2double(handles.SessionLength.String) + ...
    [0 str2double(handles.TimeWindow.String)];
set(handles.BehaviorPlot,'XLim',newLims);
set(handles.MotorPlot,'XLim',handles.SampleRate*newLims);

% Update handles structure
guidata(hObject, handles);

function TimeWindow_Callback(hObject, eventdata, handles)
Scroller_Callback(hObject, eventdata, handles);

% --- Executes when entered data in editable cell(s) in RespirationScaling.
function RespirationScaling_CellEditCallback(hObject, eventdata, handles)
handles.respiration_plot.YData = handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.SniffData';
handles.respiration_plot_gated.YData = handles.RespirationScaling.Data(1) + handles.RespirationScaling.Data(2)*handles.gatedsniffs';
guidata(hObject, handles);
