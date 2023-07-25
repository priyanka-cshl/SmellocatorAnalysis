function varargout = OEPSViewer_v4(varargin)
% OEPSVIEWER_V4 MATLAB code for OEPSViewer_v4.fig
%      OEPSVIEWER_V4, by itself, creates a new OEPSVIEWER_V4 or raises the existing
%      singleton*.
%
%      H = OEPSVIEWER_V4 returns the handle to a new OEPSVIEWER_V4 or the handle to
%      the existing singleton*.
%
%      OEPSVIEWER_V4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OEPSVIEWER_V4.M with the given input arguments.
%
%      OEPSVIEWER_V4('Property','Value',...) creates a new OEPSVIEWER_V4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OEPSViewer_v4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OEPSViewer_v4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OEPSViewer_v4

% Last Modified by GUIDE v2.5 25-Jul-2023 07:26:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OEPSViewer_v4_OpeningFcn, ...
                   'gui_OutputFcn',  @OEPSViewer_v4_OutputFcn, ...
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


% --- Executes just before OEPSViewer_v4 is made visible.
function OEPSViewer_v4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OEPSViewer_v4 (see VARARGIN)

% Choose default command line output for OEPSViewer_v4
handles.output = hObject;

% spike data 
handles.selectedUnit = plot(NaN,NaN,'w.');
handles.MySelectedUnit = patch(NaN,NaN,'w','EdgeColor','none','FaceAlpha',.8);
handles.SDline = plot(NaN,NaN,':w');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OEPSViewer_v4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OEPSViewer_v4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function WindowSize_Callback(hObject, eventdata, handles)
% hObject    handle to WindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowSize as text
%        str2double(get(hObject,'String')) returns contents of WindowSize as a double
UpdatePlot(hObject, eventdata, handles);

% --- Executes on button press in LoadSession.
function LoadSession_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.SessionPath.String)
        [RecordingFile, WhichSession] = uigetfile('/mnt/data/Sorted/Q8/2022-12-09_14-37-11/mybinaryfile.dat',...
                                'Select a processed binary file');
    handles.SessionPath.String = WhichSession;
    
    load(fullfile(handles.SessionPath.String,'chanMap.mat'),'chanMap');
    handles.NumChans.String = num2str(size(chanMap,2));
    
    % load units if available
    if exist(fullfile(WhichSession,'cluster_group.tsv')) == 2 % sssion was curated
        [handles.Units.List, handles.units.spikes] = GetSortingSummary(WhichSession);
        handles.UnitList.Data = handles.Units.List(:,[1 2 6]);
    end
end

if isempty(handles.NumChans.String)
    load(fullfile(handles.SessionPath.String,'chanMap.mat'),'chanMap');
    handles.NumChans.String = num2str(size(chanMap,2));
    
end

% get length of the file
X = dir(fullfile(handles.SessionPath.String,'mybinaryfile.dat'));
Nchan       = str2double(handles.NumChans.String);
handles.RecordingLength.String = num2str(floor(X.bytes/2/Nchan/30000/60));

% To overlay OdorTTLs
[TTLFile,TTLPath] = uigetfile('*.mat','Select a processed behavior-ephys file');
if ~isequal(TTLFile,0)
    % Load the relevant variables
    load(fullfile(TTLPath,TTLFile), 'TTLs', 'SingleUnits');
    handles.TTLs = TTLs;
    handles.SingleUnits = SingleUnits;
end

guidata(hObject, handles);

UpdatePlot(hObject, eventdata, handles);


% --- Executes on button press in ClearSession.
function UpdatePlot(hObject, eventdata, handles)
% hObject    handle to ClearSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open the binary file
OEPSSamplingRate = 30000;
Nchan       = str2double(handles.NumChans.String);
Nsamples    = str2double(handles.WindowSize.String)*OEPSSamplingRate;

TStart      = handles.TimeScroll.Value * 60 * ...
                str2double(handles.RecordingLength.String) - ...
                str2double(handles.WindowSize.String);
TStart = floor(TStart);
offset = uint64(max(0,2*Nchan * TStart * OEPSSamplingRate)); % 2 bytes per sample

fid = fopen(fullfile(handles.SessionPath.String,'mybinaryfile.dat'),'r');
fseek(fid, offset, 'bof');
MyData = fread(fid, [Nchan Nsamples], '*int16');
%MyData = fread(fid, [Nchan 33000], '*int16');
fclose(fid);

channelSpacing = str2double(handles.Spacing.String);

% subtract mean of a perticular channel from that tetrode - easy way to check for shorts
if ~isnan(handles.whichChan.Data(1,1)) && handles.chanSubtract.Value
    mychannel = handles.whichChan.Data(1,1);
    myTT = ceil(mychannel/4);
    mychannels = (myTT-1)*4 + (1:4);
    MyData(mychannels,:) = MyData(mychannels,:) - MyData(mychannel,:);
end

MyData = MyData';
MyData = fliplr(MyData);
axes(handles.axes1);
cla
hold off
MyColorOrder = brewermap(ceil(Nchan/8),'Set3')';
MyColorOrder = reshape(repmat(MyColorOrder,8,1),3,8*ceil(Nchan/8))';
set(0,'DefaultAxesColorOrder',MyColorOrder);
plot(int32(MyData) + int32(repmat(channelSpacing*(0:1:Nchan-1),size(MyData,1),1)),'LineWidth',0.11);
%set(gca,'ColorMap',brewermap(Nchan,'Accent'))
axis tight
set(gca,'Color','k','XTick',[], 'YTick',[],'YLim', [-channelSpacing Nchan*channelSpacing] );
set(0,'DefaultAxesColorOrder','remove')
myXLims = get(gca,'XLim');
myYLims = get(gca,'YLim');

myTimeAxis = max(TStart,0) + [0 str2double(handles.WindowSize.String)];
handles.myTimeAxis = myTimeAxis;

hold on
handles.selectedUnit = plot(NaN,NaN,'w.');
handles.MySelectedUnit = patch(NaN,NaN,'w','EdgeColor','none','FaceAlpha',.2);
PlotUnits_Callback(hObject, eventdata, handles);

% show SD level for selected channel
if ~isnan(handles.whichChan.Data(1,2)) && handles.chanSD.Value
    mychannel = handles.whichChan.Data(1,1);
    sd = double(channelSpacing*(size(MyData,2)-mychannel)) + (mean(MyData(:,mychannel))) - (handles.whichChan.Data(1,2)*std(double(MyData(:,mychannel))));    
    line(myXLims, [sd sd], 'LineStyle', ':', 'Color', 'w');
end

% find Trial TTLs that span this stretch
if isfield(handles,'TTLs')
    TrialStarts = find((handles.TTLs.Trial(:,1)>=myTimeAxis(1))&(handles.TTLs.Trial(:,1)<=myTimeAxis(2)));
    if ~isempty(TrialStarts)
        hold on
        Timestamps = handles.TTLs.Trial(TrialStarts,:);
        for odor = 1:3
            if any(Timestamps(:,end)==odor)
                whichTS = Timestamps(find(Timestamps(:,end)==odor),:);
                for i = 1:size(whichTS,1)
                    x(1) = whichTS(i,1)+ whichTS(i,4) - myTimeAxis(1);
                    x(2) = whichTS(i,2) - myTimeAxis(1);
                    x = x*OEPSSamplingRate; % in samples
                    Vx = [x(1) x(1) x(2) x(2)];
                    Vy = [myYLims(1) myYLims(2) myYLims(2) myYLims(1)];
                    switch odor
                        case 1
                            patch(Vx,Vy,'y','EdgeColor','none','FaceAlpha',.2);
                        case 2
                            patch(Vx,Vy,'r','EdgeColor','none','FaceAlpha',.2);
                        case 3
                            patch(Vx,Vy,'b','EdgeColor','none','FaceAlpha',.2);
                    end
                end
            end
        end
    end
end
set(gca,'XLim',myXLims);
handles.fullZoom = myYLims;
if isfield (handles, 'currentZoom')
    set(handles.axes1,'YLim',handles.currentZoom);
end
guidata(hObject, handles);

% mylims = get(gca,'XLim');
% [cluster] = GetSingleUnits(handles.SessionPath.String);
% ST = 30000*cluster(5).spikes;
% hold on
% plot(ST(:),1400,'ow');
% set(gca,'XLim',mylims);

% --- Executes on button press in ClearSession.
function ClearSession_Callback(hObject, eventdata, handles)
% hObject    handle to ClearSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SessionPath.String = '';


% --- Executes on slider movement.
function TimeScroll_Callback(hObject, eventdata, handles)
% hObject    handle to TimeScroll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
UpdatePlot(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function TimeScroll_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeScroll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function WindowSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RecordingLength_Callback(hObject, eventdata, handles)
% hObject    handle to RecordingLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RecordingLength as text
%        str2double(get(hObject,'String')) returns contents of RecordingLength as a double


% --- Executes during object creation, after setting all properties.
function RecordingLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RecordingLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Spacing_Callback(hObject, eventdata, handles)
% hObject    handle to Spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Spacing as text
%        str2double(get(hObject,'String')) returns contents of Spacing as a double
UpdatePlot(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function Spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in UnitList.
function UnitList_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to UnitList (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
% if handles.PlotUnits.Value
%     UpdatePlot(hObject, eventdata, handles);
%     % plot units if desired
handles.whichUnit.String = num2str(eventdata.Indices(1));
guidata(hObject, handles);
PlotUnits_Callback(hObject, eventdata, handles)


% end

% --- Executes on button press in PlotUnits.
function PlotUnits_Callback(hObject, eventdata, handles)
% hObject    handle to PlotUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotUnits
% if handles.PlotUnits.Value
%     UpdatePlot(hObject, eventdata, handles);
% end

myTimeAxis = handles.myTimeAxis;
OEPSSamplingRate = 30000;
Nchan       = str2double(handles.NumChans.String);
channelSpacing = str2double(handles.Spacing.String);

if handles.PlotUnits.Value
    if ~isnan(str2double(handles.whichUnit.String))
        whichUnit = str2double(handles.whichUnit.String);
        spikes = handles.units.spikes{whichUnit};
        myspikes = OEPSSamplingRate*(spikes(intersect(find(spikes>=myTimeAxis(1)),find(spikes<=myTimeAxis(2)))) - myTimeAxis(1));
        
        whichTT = floor(handles.UnitList.Data(whichUnit,2));
        whichChan = (whichTT-1)*4 + 10*(handles.UnitList.Data(whichUnit,2) - whichTT);
        offset = channelSpacing*(Nchan-whichChan-1) - channelSpacing/2;
        
        % for plotting as dots
        handles.selectedUnit.XData = myspikes;
        handles.selectedUnit.YData = 0*myspikes + offset;

        % plot as a shaded box
        myspikes(:,2:3) = myspikes + OEPSSamplingRate*[-0.001 0.001]; % 1 ms window around the spike
        myspikes(:,1) = [];
%         Vx = [myspikes(:,1) myspikes(:,1) myspikes(:,2) myspikes(:,2)];
%         Vy = repmat([-channelSpacing  channelSpacing*(Nchan-1) channelSpacing*(Nchan-1) -channelSpacing],size(myspikes,1),1);
        
        offsets =  [-channelSpacing  channelSpacing*Nchan];
        TS = myspikes';
        if ~isempty(TS)
            handles.MySelectedUnit.Vertices = [ ...
                reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
                repmat([offsets(1) offsets(2) offsets(2) offsets(1)]',size(TS,2),1)];
            handles.MySelectedUnit.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
        else
            handles.MySelectedUnit.Vertices = [];
            handles.MySelectedUnit.Faces = [];
        end
    end
end


% --- Executes when entered data in editable cell(s) in whichChan.
function whichChan_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to whichChan (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
UpdatePlot(hObject, eventdata, handles);


% --- Executes on button press in showAllUnits.
function showAllUnits_Callback(hObject, eventdata, handles)
% hObject    handle to showAllUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showAllUnits


function ChanZoom_Callback(hObject, eventdata, handles)
% hObject    handle to ChanZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ChanZoom as text
%        str2double(get(hObject,'String')) returns contents of ChanZoom as a double
channelSpacing      = str2double(handles.Spacing.String);
Nchan               = str2double(handles.NumChans.String);
chLocs              = fliplr(int32(channelSpacing*(0:1:Nchan-1)));
selChans            = eval(handles.ChanZoom.String);
handles.currentZoom = fliplr(chLocs(selChans) + int32([channelSpacing -channelSpacing]));
set(handles.axes1,'YLim',handles.currentZoom);
guidata(hObject, handles);
