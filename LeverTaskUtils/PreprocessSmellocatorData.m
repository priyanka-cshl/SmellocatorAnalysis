%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = PreprocessSmellocatorData(MyFilePath)

%% Add relevant repositories
Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));
addpath(genpath(fullfile(Paths.Code,'MatlabUtils')));

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

%% check if the preprocessed version already exists - locally or on the server
savepath = fullfile(Paths.Grid.Behavior_processed,AnimalName,[MyFileName,'_processed.mat']);
if exist(savepath)
    reply = input('This session has already been processed. \nDo you want to overwrite? Y/N [Y]: ','s');
    if strcmp(reply,'N')
        return;
    end
end

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
FileLocations.Behavior = MyFilePath;
[FilePaths, MyFileName] = fileparts(MyFilePath);
disp(MyFileName);

%% Parse into trials
[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags);
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

%% Check if passive tuning was done
MyTuningTrials = []; TrialSequence = []; PassiveReplayTraces = [];
[TuningFile] = WhereTuningFile(FilePaths,MyFileName);
if ~isempty(TuningFile)
    [MyTuningTrials, TrialSequence, PassiveReplayTraces] = ParseTuningSession(TuningFile);
    disp(['Found Tuning File: ',TuningFile]);
    FileLocations.Tuning = TuningFile;
end

%% Get info from the OEPS files if available
[myephysdir] = WhereOEPSFile(MyFileName,FilePaths); % returns empty if no recording file was found
TTLs = []; ReplayTTLs = []; TuningTTLs = [];
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

if isempty(TTLs)
    disp('no matching recording file found');
else
    FileLocations.OEPS = myephysdir;
end

%% Get spikes - label spikes by trials
SingleUnits = [];
if ~isempty(TTLs)
    foo = regexp(myephysdir,[filesep,AnimalName,filesep],'split');
    myspikesdir = fullfile(Paths.Grid.Ephys_processed,AnimalName,fileparts(foo{end}));
    if exist(myspikesdir)
        FileLocations.Spikes = myspikesdir;
        SingleUnits = GetSingleUnits(myspikesdir);
        [SingleUnits] = Spikes2Trials(SingleUnits, TTLs.Trial(1:size(TrialInfo.TrialID,2),:), TuningTTLs);    
    end
end

%% Saving stuff in one place
save(savepath, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', 'FileLocations', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');
    
end