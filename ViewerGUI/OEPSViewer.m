function varargout = OEPSViewer(varargin)
% OEPSVIEWER MATLAB code for OEPSViewer.fig
%      OEPSVIEWER, by itself, creates a new OEPSVIEWER or raises the existing
%      singleton*.
%
%      H = OEPSVIEWER returns the handle to a new OEPSVIEWER or the handle to
%      the existing singleton*.
%
%      OEPSVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OEPSVIEWER.M with the given input arguments.
%
%      OEPSVIEWER('Property','Value',...) creates a new OEPSVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OEPSViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OEPSViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OEPSViewer

% Last Modified by GUIDE v2.5 18-Aug-2022 22:31:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OEPSViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @OEPSViewer_OutputFcn, ...
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


% --- Executes just before OEPSViewer is made visible.
function OEPSViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OEPSViewer (see VARARGIN)

% Choose default command line output for OEPSViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OEPSViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OEPSViewer_OutputFcn(hObject, eventdata, handles) 
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
        [RecordingFile, WhichSession] = uigetfile('*.dat',...
                                'Select a processed binary file');
    handles.SessionPath.String = WhichSession;
    
    load(fullfile(handles.SessionPath.String,'chanMap.mat'),'chanMap');
    handles.NumChans.String = num2str(size(chanMap,2));
end

if isempty(handles.NumChans.String)
    load(fullfile(handles.SessionPath.String,'chanMap.mat'),'chanMap');
    handles.NumChans.String = num2str(size(chanMap,2));
    
end

% get length of the file
X = dir(fullfile(handles.SessionPath.String,'mybinaryfile.dat'));
Nchan       = str2double(handles.NumChans.String);
handles.RecordingLength.String = num2str(floor(X.bytes/2/Nchan/32000/60));

UpdatePlot(hObject, eventdata, handles);


% --- Executes on button press in ClearSession.
function UpdatePlot(hObject, eventdata, handles)
% hObject    handle to ClearSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open the binary file
Nchan       = str2double(handles.NumChans.String);
Nsamples    = str2double(handles.WindowSize.String)*32000;

TStart      = handles.TimeScroll.Value * 60 * ...
                str2double(handles.RecordingLength.String) - ...
                str2double(handles.WindowSize.String);
offset = max(0,Nchan * TStart * 32000);

fid = fopen(fullfile(handles.SessionPath.String,'mybinaryfile.dat'),'r');
fseek(fid, offset, 'bof');
MyData = fread(fid, [Nchan 33000], '*int16');
fclose(fid);

channelSpacing = str2double(handles.Spacing.String);

MyData = MyData';
axes(handles.axes1);
MyColorOrder = brewermap(ceil(Nchan/8),'Set3')';
MyColorOrder = reshape(repmat(MyColorOrder,8,1),3,8*ceil(Nchan/8))';
set(0,'DefaultAxesColorOrder',MyColorOrder);
plot(MyData + int16(repmat(channelSpacing*(0:1:Nchan-1),size(MyData,1),1)));
%set(gca,'ColorMap',brewermap(Nchan,'Accent'))
axis tight
set(gca,'Color','k','XTick',[], 'YTick',[]);
set(0,'DefaultAxesColorOrder','remove')
guidata(hObject, handles);

% --- Executes on button press in ClearSession.
function ClearSession_Callback(hObject, eventdata, handles)
% hObject    handle to ClearSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
