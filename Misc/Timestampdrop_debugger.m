%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = Timestampdrop_debugger(MyFilePath,overwriteflag)

%% Add relevant repositories
Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));
addpath(genpath(fullfile(Paths.Code,'open-ephys-matlab-tools'))); % for the new OEPS GUI: https://github.com/open-ephys/open-ephys-matlab-tools (use commit 10-04-22)
addpath(genpath(fullfile(Paths.Code,'MatlabUtils')));
addpath(genpath(fullfile(Paths.Code,'npy-matlab/'))); % path to npy-matlab scripts

%% globals
% global MyFileName;
global SampleRate;
SampleRate = 500; % Samples/second
global startoffset;
startoffset = 1; % in seconds
global savereplayfigs;
savereplayfigs = 0;
global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips]
errorflags = [0 0 0 0];
global TargetZones; %#ok<*NUSED>

%% core data extraction (and settings)
if ~exist(MyFilePath)
    foo = regexp(MyFilePath,'_','split');
    AnimalName = foo{1};
    MyFilePath = fullfile(Paths.Grid.Behavior,AnimalName,MyFilePath);
    [FilePaths, MyFileName] = fileparts(MyFilePath); %#ok<*ASGLU>
else
    [FilePaths, MyFileName] = fileparts(MyFilePath);
    [~,AnimalName] = fileparts(FilePaths);
end

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);

% read data into chunks as they came in
global samplenum; % #samples read in the current session
global TotalData; % matrix containing data from the current callback (refreshed at 20Hz (or Refresh rate))
global TotalTime; % matrix containing timestamps from the current callback
samplenum = 0; TotalData = zeros(6000,1); TotalTime = zeros(6000,1);
samps2read = 25;

handles.axes1 = axes; % main plot
handles.trial_on_1 = fill(NaN,NaN,[.8 .8 .8]);
hold on;
handles.trial_on_1.EdgeColor = 'none';
handles.trial_on_2 = fill(NaN,NaN,[0.8941    0.9412    0.9020]);
handles.trial_on_2.EdgeColor = 'none';
handles.trial_on_3 = fill(NaN,NaN,[0.8706    0.9216    0.9804]);
handles.trial_on_3.EdgeColor = 'none';
handles.trial_on_4 = fill(NaN,NaN,[0.93    0.84    0.84]);
handles.trial_on_4.EdgeColor = 'none';

while (samplenum + 25) <= size(MyData,1)
    newsamps = MyData(samplenum+(1:samps2read),:);
    last_ts = TotalTime(end,1);
    first_ts = newsamps(1,1);
    if first_ts - last_ts > 0.003
        keyboard;
        newsamps(:,1) = last_ts + (0.002:0.002:0.05);
        newsamps(:,6) = [zeros(24,1); MyData(samplenum+26,6)];
    else
        samplenum = samplenum + 25;
    end
    TotalTime = [ TotalTime(samps2read+1:end,1); newsamps(:,1) ];
    TotalData = [ TotalData(samps2read+1:end,1); newsamps(:,6) ];
    [handles] = PlotToPatch_Trial_GUI(handles, TotalData(:,1), TotalTime, [0 5]);

end


%% Parse into trials
%[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[Trials,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

% sanity check - did some guess work in CorrectMatlabSampleDrops to compute
% odor start - check if it made sense
if ~isempty(InitiationsFixed)
    if any(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01)
        weirdo = find(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01);
        if any(TrialInfo.OdorStart(InitiationsFixed(weirdo),2)>-1)
            disp('something funky with computing odorstart from Lever trace');
            keyboard;
            TrialInfo.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo.OdorStart(InitiationsFixed(weirdo),2);
        else
            % Initiation hold was larger than a second - that's couldn't
            % compute it accurately from trial traces
            TrialInfo.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo.OdorStart(InitiationsFixed(weirdo),2);
        end
    end
end

%% Check if passive tuning was done
MyTuningTrials = []; TrialSequence = []; PassiveReplayTraces = [];
[TuningFile] = WhereTuningFile(FilePaths,MyFileName);
if ~isempty(TuningFile)
    [MyTuningTrials, TrialSequence, PassiveReplayTraces] = ParseTuningSession(TuningFile);
    disp(['Found Tuning File: ',TuningFile]);
    FileLocations.Tuning = TuningFile;
end

%% Get info from the OEPS files if available
TTLs = []; ReplayTTLs = []; TuningTTLs = []; myephysdir = [];

% for batch Q - first check if event data has been processed during Sorting
[mySortingdir] = WhereOEPSTTLs(MyFileName,FilePaths); % returns empty if no recording file was found
if ~isempty(mySortingdir)
    % process the TTLs
    [TTLs,ReplayTTLs,TuningTTLs,~] = ...
            GetOepsAuxChannels(mySortingdir, Trials.TimeStamps, MyTuningTrials, TrialSequence, 'KiloSorted', 1); % send 'ADC', 1 to also get analog aux data
end

if isempty(TTLs)
    [myephysdir] = WhereOEPSFile(MyFileName,FilePaths); % returns empty if no recording file was found
    if ~isempty(myephysdir)
        if size(myephysdir,1) == 1
            [TTLs,ReplayTTLs,TuningTTLs,~] = ...
                GetOepsAuxChannels(myephysdir, Trials.TimeStamps, MyTuningTrials, TrialSequence); % send 'ADC', 1 to also get analog aux data
        else
            while isempty(TTLs) && ~isempty(myephysdir)
                [TTLs,ReplayTTLs,TuningTTLs,~] = ...
                    GetOepsAuxChannels(myephysdir(1,:), Trials.TimeStamps, MyTuningTrials, TrialSequence);
                if isempty(TTLs)
                    myephysdir(1,:) = [];
                end
            end
        end
    end
else
    FileLocations.KiloSorted = mySortingdir;
end

if isempty(TTLs)
    disp('no matching recording file found');
else
    FileLocations.OEPS = myephysdir;
end

%% Get spikes - label spikes by trials
SingleUnits = [];
if ~isempty(TTLs)
    if isfield(FileLocations,'KiloSorted')
        myspikesdir = FileLocations.KiloSorted;
    else
        foo = regexp(myephysdir,[filesep,AnimalName,filesep],'split');
        myspikesdir = fullfile(Paths.Local.Ephys_processed,AnimalName,fileparts(foo{end}));
    end
    if exist(myspikesdir)
        FileLocations.Spikes = myspikesdir;
        SingleUnits = GetSingleUnits(myspikesdir);
        [SingleUnits] = Spikes2Trials(SingleUnits, TTLs.Trial(1:size(TrialInfo.TrialID,2),:), TuningTTLs);    
    end
end

%% Photometry
FTrace = [];
if ~isempty(TTLs)
    
end

%% Saving stuff in one place
if ~exist(fileparts(savepath),'dir')
    mkdir(fileparts(savepath));
end
save(savepath, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', 'FileLocations', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');
    
end