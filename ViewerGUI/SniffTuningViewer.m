function varargout = SniffTuningViewer(varargin)
% SNIFFTUNINGVIEWER MATLAB code for SniffTuningViewer.fig
%      SNIFFTUNINGVIEWER, by itself, creates a new SNIFFTUNINGVIEWER or raises the existing
%      singleton*.
%
%      H = SNIFFTUNINGVIEWER returns the handle to a new SNIFFTUNINGVIEWER or the handle to
%      the existing singleton*.
%
%      SNIFFTUNINGVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNIFFTUNINGVIEWER.M with the given input arguments.
%
%      SNIFFTUNINGVIEWER('Property','Value',...) creates a new SNIFFTUNINGVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SniffTuningViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SniffTuningViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SniffTuningViewer

% Last Modified by GUIDE v2.5 09-Jul-2025 11:18:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SniffTuningViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @SniffTuningViewer_OutputFcn, ...
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


% --- Executes just before SniffTuningViewer is made visible.
function SniffTuningViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SniffTuningViewer (see VARARGIN)

% Choose default command line output for SniffTuningViewer
handles.output = hObject;

handles.figure1.Position(4) = 35;

%% make TABs
%Create tab groupA - sniff locked vs. contamination
handles.tgroupA = uitabgroup('Parent', handles.figure1,'TabLocation', 'top',...
    'Position', [0.145 0.09 0.8500 0.8500] );
handles.tabA1 = uitab('Parent', handles.tgroupA, 'Title', 'Stimulus-wise');
handles.tabA2 = uitab('Parent', handles.tgroupA, 'Title', 'Contanimation');
%Place panels into each tab
set(handles.A1,'Parent',handles.tabA1)
set(handles.A2,'Parent',handles.tabA2)
%Reposition each panel to same location as panel 1
%handles.A1.Position = [1.000    0.1   172.000   28.0000];
handles.A1.Position = [0.01    0.01   0.98   0.97];
set(handles.A2,'position',get(handles.A1,'position'));

%Create tab groupB - single unit or multi unit
handles.tgroupB = uitabgroup('Parent', handles.figure1,'TabLocation', 'top',...
    'Position', [0.005 0.38 0.135 0.560] );
handles.tabB1 = uitab('Parent', handles.tgroupB, 'Title', 'Basic');
handles.tabB2 = uitab('Parent', handles.tgroupB, 'Title', 'Multi');
%Place panels into each tab
set(handles.View1,'Parent',handles.tabB1)
set(handles.View2,'Parent',handles.tabB2)
%Reposition each panel to same location as panel 1
handles.View1.Position = [0.05    0.1   27.000   17.0000];
set(handles.View2,'position',get(handles.View1,'position'));

%% session path
[Paths] = WhichComputer();
handles.quickmode = 0;
handles.AllUnits = [];
handles.RefreshUnits.Enable = 'off';
handles.snifftypes = 5;

if ~isempty(varargin)
    if numel(varargin) == 2 && strcmp(varargin{2},'quickprocess')
        if exist(varargin{1}) == 7
            handles.WhereSession.String = varargin{1};   
        else
            handles.WhereSession.String = '/mnt/data/Sorted/NS12/2023-07-27_17-18-37';
        end
        [~,handles.MouseName] = fileparts(fileparts(handles.WhereSession.String));
        handles.quickmode = 1;
        handles.RefreshUnits.Enable = 'on';
    else
        if exist(varargin{1}) == 2
            handles.WhereSession.String = varargin{1};
        else
            handles.WhereSession.String = fullfile(Paths.Wolf.processed,'forWDW',varargin{1});
        end
        [~,FileName] = fileparts(handles.WhereSession.String);
        handles.MouseName = regexprep(FileName,'_(\w+)_processed.mat','');
    end
else
    handles.WhereSession.String = fullfile(Paths.Wolf.processed,'forWDW','Q4_20221112_r0_processed.mat');
    [~,FileName] = fileparts(handles.WhereSession.String);
    handles.MouseName = regexprep(FileName,'_(\w+)_processed.mat','');
end


%% initializations
for i = 1:5
    currentaxis = handles.(['axes',num2str(i)]);
    axes(currentaxis);
    set(currentaxis,'Color','none','XTick',[],'YTick',[]);
    hold on
    handles.(['overlayspikesplot',num2str(i)]) = plot(NaN, NaN, '.','Color', [0.7 0.5 0.85], 'Markersize', 0.75);
    handles.(['spikesplot',num2str(i)]) = plot(NaN, NaN, '.k','Markersize', 0.5);
end
handles.stackMerge.Enable = 'off';

% Update handles structure
guidata(hObject, handles);

if exist(handles.WhereSession.String)==2 || exist(handles.WhereSession.String)==7
    handles = LoadSession_Callback(hObject, eventdata, handles);
end

handles.figure1.Resize = 'on';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SniffTuningViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SniffTuningViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadSession.
function [handles] = LoadSession_Callback(hObject, eventdata, handles)
    
    if ~handles.quickmode
        if isempty(handles.WhereSession.String)
            [Paths] = WhichComputer();
            [WhichSession, SessionPath] = uigetfile(...
                fullfile(Paths.Wolf.processed,'forWDW','Q4_20221112_r0_processed.mat'),...
                'Select Session');
            handles.WhereSession.String = fullfile(SessionPath,WhichSession);
        end

        % get a list of sniffs, also get units
        [handles.AllSniffs, handles.SniffColumnInfo, SingleUnits] = GetAllSniffs(handles.WhereSession.String);
    else
        % get a list of sniffs, also get units
        [handles.AllSniffs] = QuickSniffTTLMapper_v2(handles.WhereSession.String);
        if numel(unique(handles.AllSniffs(:,4))) == 1
            handles.snifftypes = 1;
            handles.stack_many_units.Value = 0;
        end

        % loading units
        SingleUnits = GetSingleUnits(handles.WhereSession.String, 3);
    end
    
    handles.MultiSelectUnits = [];
    handles.selectCount.Data(1) = 0;
    for whichunit = 1:size(SingleUnits,2)
        handles.AllUnits.Spikes{whichunit}      = SingleUnits(whichunit).spikes; % raw timestamps in OEPS base
        handles.AllUnits.ChannelInfo(whichunit,1:2) = [SingleUnits(whichunit).tetrode SingleUnits(whichunit).id]; % tetrode and phy cluster ID
        handles.UnitList.Data(whichunit,:) = [whichunit handles.AllUnits.ChannelInfo(whichunit,[2 1])];
    end
    
    handles.NumUnits.Data(1) = size(handles.AllUnits.Spikes,2);
    if isnan(handles.thisUnit.Data(1,1)) || handles.thisUnit.Data(1,1)>handles.NumUnits.Data(1)
        handles.thisUnit.Data(1,1) = 1;
    end

    handles = UpdatePlots(hObject, eventdata, handles);

    % Update handles structure
    guidata(hObject, handles);


% --- Executes within Load Session and at any unit update call.
function [handles] = UpdatePlots(hObject, eventdata, handles)

    % get sniff time stamps and info for the sniffs we want to plot
    [handles.SelectedSniffs] = ParseSniffsByType(handles.AllSniffs, [], 'SortBy', handles.SortSniffsBy.Value);
    %myXlim = eval(handles.xlims.String);
    
    handles.maxsniffs = max(cellfun(@length, handles.SelectedSniffs)) + 100*(max(handles.SelectedSniffs{1}(:,8))-1);
    for snifftype = 1:handles.snifftypes
        set(handles.(['axes',num2str(snifftype)]), 'YLim',[0 handles.maxsniffs+1]);
    end

    % Update handles structure
    guidata(hObject, handles);

    UpdateUnits(handles);
    if handles.CompareUnits.Value
        CompareUnits_Callback(hObject, eventdata, handles);
    end

% --- Executes within Load Session and at any unit update or plot update call.
function UpdateUnits(handles)

    whichUnit = handles.thisUnit.Data(1,1);
    handles.thisUnit.Data(2,1) = handles.AllUnits.ChannelInfo(whichUnit,2); % cluster ID
    handles.thisUnit.Data(3,1) = handles.AllUnits.ChannelInfo(whichUnit,1); % tetrode
    %myXlim = eval(handles.xlims.String);
    
    thisUnitSpikes = handles.AllUnits.Spikes{whichUnit};
    
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(handles.SelectedSniffs, thisUnitSpikes);
    if handles.PoolUnits.Value 
        PoolUnitSpikes = handles.AllUnits.Spikes{handles.thisUnit.Data(1,2)};
        [Pool2Raster, ~] = GetSniffLockedSpikes(handles.SelectedSniffs, PoolUnitSpikes);
        for snifftype = 1:handles.snifftypes
            SpikeRaster{snifftype} = vertcat(SpikeRaster{snifftype},Pool2Raster{snifftype});
        end
    end
    
    for snifftype = 1:handles.snifftypes
        myspikeplot = ['spikesplot',num2str(snifftype)];
        SpikesPlot = SpikeRaster{snifftype};
        if ~isempty(SpikesPlot)
            SpikesPlot(:,2) = abs(SpikesPlot(:,2));
            set(handles.(myspikeplot),'XData',SpikesPlot(:,1),'YData',SpikesPlot(:,2));
        end
    end


% --- Executes on button press in NextUnit.
function NextUnit_Callback(hObject, eventdata, handles)
    temp = rem(handles.thisUnit.Data(1,1)+1, handles.NumUnits.Data(1));
    if ~temp
        handles.thisUnit.Data(1,1) = handles.NumUnits.Data(1);
    else
        handles.thisUnit.Data(1,1) = temp;
    end
    UpdateUnits(handles);
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in PrevUnit.
function PrevUnit_Callback(hObject, eventdata, handles)
    temp = rem(handles.thisUnit.Data(1,1)-1, handles.NumUnits.Data(1));
    if ~temp
        handles.thisUnit.Data(1,1) = handles.NumUnits.Data(1);
    else
        handles.thisUnit.Data(1,1) = temp;
    end
    UpdateUnits(handles);
    % Update handles structure
    guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in thisUnit.
function thisUnit_CellEditCallback(hObject, eventdata, handles)
    if eventdata.Indices(1) == 2
        col = eventdata.Indices(2);
        if ~isempty(find(handles.AllUnits.ChannelInfo(:,2)==eventdata.NewData))
            handles.thisUnit.Data(1,col) = find(handles.AllUnits.ChannelInfo(:,2)==eventdata.NewData);
            handles.thisUnit.Data(3,col) = handles.AllUnits.ChannelInfo(handles.thisUnit.Data(1,col),1);
        end
    end

    if eventdata.Indices(1) == 1 && eventdata.Indices(2) == 2
        col = eventdata.Indices(2);
        if ~isempty(find(handles.AllUnits.ChannelInfo(:,2)==eventdata.NewData))
            handles.thisUnit.Data(2,col) = handles.AllUnits.ChannelInfo(handles.thisUnit.Data(1,col),2);
            handles.thisUnit.Data(3,col) = handles.AllUnits.ChannelInfo(handles.thisUnit.Data(1,col),1);
        end
    end

    % Update handles structure
    guidata(hObject, handles);
    UpdateUnits(handles);
    % Update handles structure
    guidata(hObject, handles);

% --- Executes when entered data in editable cell(s) in ComparedUnit.
function ComparedUnit_CellEditCallback(hObject, eventdata, handles)
    if eventdata.Indices(1) == 2
        if ~isempty(find(handles.AllUnits.ChannelInfo(:,2)==eventdata.NewData))
            handles.ComparedUnit.Data(1) = find(handles.AllUnits.ChannelInfo(:,2)==eventdata.NewData);
            handles.ComparedUnit.Data(3) = handles.AllUnits.ChannelInfo(handles.ComparedUnit.Data(1),1);
        end
    end
    % Update handles structure
    guidata(hObject, handles);
    if handles.CompareUnits.Value
        CompareUnits_Callback(hObject, eventdata, handles)
    end


% --- Executes on selection change in SortSniffsBy.
function SortSniffsBy_Callback(hObject, eventdata, handles)
    handles = UpdatePlots(hObject, eventdata, handles);
    % Update handles structure
    guidata(hObject, handles);

    % --- Executes on button press in PoolUnits.
function PoolUnits_Callback(hObject, eventdata, handles)
    UpdateUnits(handles);

% --- Executes on button press in CompareUnits.
function CompareUnits_Callback(hObject, eventdata, handles)
    if ~handles.CompareUnits.Value
        handles.stackMerge.Enable = 'off';
        for snifftype = 1:handles.snifftypes
            set(handles.(['overlayspikesplot',num2str(snifftype)]),'XData',nan,'YData',nan);
        end
    else
        handles.stackMerge.Enable = 'on';
        OverlayUnit = handles.ComparedUnit.Data(1);
        handles.ComparedUnit.Data(3) = handles.AllUnits.ChannelInfo(OverlayUnit,1);
        handles.ComparedUnit.Data(2) = handles.AllUnits.ChannelInfo(OverlayUnit,2);

        MergeUnitSpikes = handles.AllUnits.Spikes{OverlayUnit};
        [OverlaySpikeRaster, ~] = GetSniffLockedSpikes(handles.SelectedSniffs, MergeUnitSpikes);

        for snifftype = 1:handles.snifftypes
            myspikeplot = ['overlayspikesplot',num2str(snifftype)];
            SpikesPlot = OverlaySpikeRaster{snifftype};
            SpikesPlot(:,2) = abs(SpikesPlot(:,2));
%             if handles.stackMerge.Value
%                 set(handles.(myspikeplot),'XData',SpikesPlot(:,1),'YData', ...
%                     SpikesPlot(:,2) + handles.maxsniffs + 100 );
%             else
                set(handles.(myspikeplot),'XData',SpikesPlot(:,1),'YData',SpikesPlot(:,2));
%             end
        end
    end
    stackMerge_Callback(hObject, eventdata, handles);
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in stackMerge.
function stackMerge_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of stackMerge
    if handles.stackMerge.Value && handles.CompareUnits.Value % stack the two merging units
        % change ylim
        ylim = [-(handles.maxsniffs + 1) (handles.maxsniffs + 1)];
        spikecolor  = [0.3 0.05 0.2]; %[0 0.5 0.75];
        spikesize = 0.5;
%         handles.figure1.Position(4) = 35*2;
%         handles.A1.Position(4) = 28*2;
    else
        % change ylim
        ylim = [0 (handles.maxsniffs + 1)];
        spikecolor = [0.7000 0.5000 0.8500];
        spikesize = 0.75;
%         handles.figure1.Position(4) = 35;
%         handles.A1.Position(4) = 28;
    end
    for i = 1:5
        set(handles.(['axes',num2str(i)]),'YLim', ylim);
        myspikeplot = ['overlayspikesplot',num2str(i)];
        myYData = abs(get(handles.(myspikeplot),'YData'));
        if handles.stackMerge.Value
            myYData = -myYData;
        end
        set(handles.(myspikeplot),'YData',myYData,'Color',spikecolor,'Markersize', spikesize);
    end
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in RefreshUnits.
function RefreshUnits_Callback(hObject, eventdata, handles)
    % loading units
    SingleUnits = GetSingleUnits(handles.WhereSession.String, 3);
    handles.AllUnits = [];
    handles.MultiSelectUnits = [];
    handles.selectCount.Data(1) = 0;
    for whichunit = 1:size(SingleUnits,2)
        handles.AllUnits.Spikes{whichunit}      = SingleUnits(whichunit).spikes; % raw timestamps in OEPS base
        handles.AllUnits.ChannelInfo(whichunit,1:2) = [SingleUnits(whichunit).tetrode SingleUnits(whichunit).id]; % tetrode and phy cluster ID
        handles.UnitList.Data(whichunit,:) = [whichunit handles.AllUnits.ChannelInfo(whichunit,[2 1])];
    end
    if size(handles.UnitList.Data,1) > whichunit
        handles.UnitList.Data(whichunit:end,:) = [];
    end
    
    handles.NumUnits.Data(1) = size(handles.AllUnits.Spikes,2);
    if isnan(handles.thisUnit.Data(1,1)) || handles.thisUnit.Data(1,1)>handles.NumUnits.Data(1)
        handles.thisUnit.Data(1,1) = 1;
    end

    % Update handles structure
    guidata(hObject, handles);
    
    UpdateUnits(handles);


% --- Executes on button press in ShowMulti.
function ShowMulti_Callback(hObject, eventdata, handles)
    ygap = 500;
    if handles.selectCount.Data(1)
        for n = 1:handles.selectCount.Data(1)
            if handles.snifftypes == 1 && ~handles.stack_many_units.Value
                whichUnit = handles.MultiSelectUnits(n);
            else
                whichUnit = handles.MultiSelectUnits(handles.selectCount.Data(1) + 1 - n);
            end
            thisUnitSpikes = handles.AllUnits.Spikes{whichUnit};
            if n == 1
                [PooledRaster, maxsniffs] = GetSniffLockedSpikes(handles.SelectedSniffs, thisUnitSpikes);
                % maxsniffs = maxsniffs + 100;
                maxsniffs = handles.maxsniffs;
            else
                if handles.snifftypes == 1 && ~handles.stack_many_units.Value
                    if n <= 5
                        [SpikeRaster, ~] = GetSniffLockedSpikes(handles.SelectedSniffs, thisUnitSpikes);
                        PooledRaster{n} = SpikeRaster{1};
                        stackjump = 1;
                    else
                        whichplot = rem(n-1,5)+1;
                        stackjump = floor((n-1)/5) + 1;
                        [SpikeRaster, ~] = GetSniffLockedSpikes(handles.SelectedSniffs, thisUnitSpikes, ((stackjump-1)*(maxsniffs + ygap)));
                        PooledRaster{whichplot} = vertcat(SpikeRaster{1},PooledRaster{whichplot});
                    end
                else
                    [SpikeRaster, ~] = GetSniffLockedSpikes(handles.SelectedSniffs, thisUnitSpikes, ((n-1)*(maxsniffs + ygap)));
                    for snifftype = 1:handles.snifftypes
                        PooledRaster{snifftype} = vertcat(SpikeRaster{snifftype},PooledRaster{snifftype});
                    end
                end
            end
        end
        
        if handles.snifftypes == 1 && ~handles.stack_many_units.Value
            for snifftype = 1:max(5,(rem(n-1,5)+1))
                myspikeplot = ['spikesplot',num2str(snifftype)];
                SpikesPlot = PooledRaster{snifftype};
                if ~isempty(SpikesPlot)
                    SpikesPlot(:,2) = abs(SpikesPlot(:,2));
                    set(handles.(myspikeplot),'XData',SpikesPlot(:,1),'YData',SpikesPlot(:,2));
                    set(handles.(['axes',num2str(snifftype)]), 'YLim',[0 stackjump*(handles.maxsniffs+ygap)]);
                end
            end

        else
            for snifftype = 1:handles.snifftypes
                myspikeplot = ['spikesplot',num2str(snifftype)];
                SpikesPlot = PooledRaster{snifftype};
                if ~isempty(SpikesPlot)
                    SpikesPlot(:,2) = abs(SpikesPlot(:,2));
                    set(handles.(myspikeplot),'XData',SpikesPlot(:,1),'YData',SpikesPlot(:,2));
                    set(handles.(['axes',num2str(snifftype)]), 'YLim',[0 n*(handles.maxsniffs+ygap)]);
                end
            end
        end
    end

% --- Executes when selected cell(s) is changed in UnitList.
function UnitList_CellSelectionCallback(hObject, eventdata, handles)
    if ~isempty(eventdata.Indices)
        selectedUnits = eventdata.Indices(:,1);
        handles.MultiSelectUnits = selectedUnits;
        handles.selectCount.Data(1) = numel(selectedUnits);
    end
    % Update handles structure
    guidata(hObject, handles);
    


% --- Executes on button press in RescaleYLim.
function RescaleYLim_Callback(hObject, eventdata, handles)
% hObject    handle to RescaleYLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    for snifftype = 1:handles.snifftypes
        set(handles.(['axes',num2str(snifftype)]), 'YLim',[0 handles.maxsniffs+1]);
    end
