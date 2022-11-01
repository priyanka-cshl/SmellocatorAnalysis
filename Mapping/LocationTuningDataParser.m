function [PassiveTuning] = LocationTuningDataParser(MyFilePath)

%% Add relevant repositories
Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));
addpath(genpath(fullfile(Paths.Code,'MatlabUtils')));

%% globals
global SampleRate;
SampleRate = 500; % Samples/second
global startoffset;
startoffset = 1; % in seconds

%% core data extraction (and settings)
if  exist(MyFilePath) == 2
    foo = regexp(MyFilePath,'_','split');
    AnimalName = foo{1};
    MyFilePath = fullfile(Paths.Mapping.StimulusFile,AnimalName,MyFilePath);
    [FilePaths, MyFileName] = fileparts(MyFilePath); %#ok<*ASGLU>
else
    [FilePaths, MyFileName] = fileparts(MyFilePath);
    [~,AnimalName] = fileparts(FilePaths);
end

%% Get Mapping Trial Sequence from the 'behavior' file
MyTuningTrials = []; TrialSequence = [];
TuningFile = MyFilePath;
if ~isempty(TuningFile)
    [MyTuningTrials, TrialSequence, ~, TuningParams] = ParseTuningSession(TuningFile);
    disp(['Found Tuning File: ',TuningFile]);
    FileLocations.Tuning = TuningFile;
end

% TuningParams: [0 121 MotorSettle PreOdor Odor Purge PostOdor ITI LeverGain LeverDC Reward]

%% Get info from the OEPS files if available
[myephysdir] = WhereOEPSTuningFile(MyFileName,FilePaths); % returns empty if no recording file was found
TTLs = []; TuningTTLs = [];
if ~isempty(myephysdir)
    if size(myephysdir,1) == 1
        [TTLs,~,TuningTTLs,~] = ...
            GetOepsAuxChannels(myephysdir, [], MyTuningTrials, TrialSequence); % send 'ADC', 1 to also get analog aux data
    else
        while isempty(TTLs) && ~isempty(myephysdir)
            [TTLs,~,TuningTTLs,~] = ...
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

%% Get all units and rasters
[TuningSpikes, TuningTTLs] = TrialAlignedSpikeTimesPassiveTuning(SingleUnits,TTLs,TuningTTLs);
window = [-TuningParams(1,8)/2 (sum(TuningParams(4:7)) + TuningParams(1,8)/2)]/1000;
PassiveTuning.FileLocations = FileLocations;
PassiveTuning.ITI = 1000*[-window(1) TuningParams(1,8)/1000+window(1)];
PassiveTuning.PreOdor = TuningParams(1,4);
PassiveTuning.Odor = TuningParams(1,5);
PassiveTuning.PostOdor = TuningParams(1,6) + TuningParams(1,7);
PassiveTuning.RasterDimensions = {'Units'; 'Odors'; 'Locations'; 'Time'; 'Repeats'};
[PassiveTuning.RasterOut,PassiveTuning.SpikeCounts,PassiveTuning.Odors,PassiveTuning.Locations] = ForPassiveTuningDecoder(TuningSpikes, TuningTTLs, window);
    
end