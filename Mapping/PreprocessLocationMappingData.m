%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = PreprocessLocationMappingData(MyFilePath)

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

%% core data extraction (and settings)
if ~exist(MyFilePath)
    foo = regexp(MyFilePath,'_','split');
    AnimalName = foo{1};
    MyFilePath = fullfile(Paths.Mapping.StimulusFile,AnimalName,MyFilePath);
    [FilePaths, MyFileName] = fileparts(MyFilePath); %#ok<*ASGLU>
else
    [FilePaths, MyFileName] = fileparts(MyFilePath);
    [~,AnimalName] = fileparts(FilePaths);
end

%% check if the preprocessed version already exists - locally or on the server
savepath = fullfile(Paths.Mapping.Processed,AnimalName,[MyFileName,'_processed.mat']);
if exist(savepath)
    reply = input('This session has already been processed. \nDo you want to overwrite? Y/N [Y]: ','s');
    if strcmp(reply,'N')
        return;
    end
end

%% Get Mapping Trial Sequence from the 'behavior' file
MyTuningTrials = []; TrialSequence = [];
TuningFile = MyFilePath;
if ~isempty(TuningFile)
    [MyTuningTrials, TrialSequence] = ParseTuningSession(TuningFile);
    disp(['Found Tuning File: ',TuningFile]);
    FileLocations.Tuning = TuningFile;
end

%% Get info from the OEPS files if available
[myephysdir] = WhereOEPSTuningFile(MyFileName,FilePaths); % returns empty if no recording file was found
TTLs = []; ReplayTTLs = []; TuningTTLs = [];
if ~isempty(myephysdir)
    if size(myephysdir,1) == 1
        [TTLs,ReplayTTLs,TuningTTLs,~] = ...
            GetOepsAuxChannels(myephysdir, [], MyTuningTrials, TrialSequence); % send 'ADC', 1 to also get analog aux data
    else
        while isempty(TTLs) && ~isempty(myephysdir)
            [TTLs,ReplayTTLs,TuningTTLs,~] = ...
                GetOepsAuxChannels(myephysdir(1,:), [], MyTuningTrials, TrialSequence);
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
    myspikesdir = fullfile(Paths.Mapping.EphysSorted,AnimalName,fileparts(foo{end}));
    if exist(myspikesdir)
        FileLocations.Spikes = myspikesdir;
        SingleUnits = GetSingleUnits(myspikesdir);
        [SingleUnits] = Spikes2Trials(SingleUnits, [], TuningTTLs);    
    end
end

%% Saving stuff in one place
save(savepath, 'startoffset', 'SampleRate', 'FileLocations', ...
               'TTLs', 'TuningTTLs', 'SingleUnits');
    
end