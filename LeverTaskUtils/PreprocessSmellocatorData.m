%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = PreprocessSmellocatorData(MyFilePath,overwriteflag)

if nargin<2
    overwriteflag = 0;
end

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
if exist(MyFilePath) & ~isempty(fileparts(MyFilePath))
    [FilePaths, MyFileName] = fileparts(MyFilePath);
    [~,AnimalName] = fileparts(FilePaths);
else
    foo = regexp(MyFilePath,'_','split');
    AnimalName = foo{1};
    MyFilePath = fullfile(Paths.Grid.Behavior,AnimalName,MyFilePath);
    [FilePaths, MyFileName] = fileparts(MyFilePath); %#ok<*ASGLU>
end

%% check if the preprocessed version already exists - locally or on the server
savepath = fullfile(Paths.Grid.Behavior_processed,AnimalName,[MyFileName,'_processed.mat']);
if ~overwriteflag && exist(savepath)
    reply = input('This session has already been processed. \nDo you want to overwrite? Y/N [Y]: ','s');
    if strcmp(reply,'N')
        return;
    end
end

%[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath,'Fastmode',1);
FileLocations.Behavior = MyFilePath;
[FilePaths, MyFileName] = fileparts(MyFilePath);
disp(MyFileName);

%% Parse into trials
%[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[Trials,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

% sanity check - did some guess work in CorrectMatlabSampleDrops to compute
% odor start - check if it made sense
problematic_odor_starts = [];
if ~isempty(InitiationsFixed)
    if any(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01)
        weirdo = find(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01);
        problematic_odor_starts = InitiationsFixed(weirdo);
%         if any(TrialInfo.OdorStart(InitiationsFixed(weirdo),2)>-1)
%             disp('something funky with computing odorstart from Lever trace');
%             keyboard;
%             if ~errorflags(1)
%                 TrialInfo.OdorStart(InitiationsFixed(weirdo),2) = TrialInfo.OdorStart(InitiationsFixed(weirdo),1);
%             else
%                 TrialInfo.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo.OdorStart(InitiationsFixed(weirdo),2);
%             end
%         else
%             % Initiation hold was larger than a second - that's couldn't
%             % compute it accurately from trial traces
%             TrialInfo.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo.OdorStart(InitiationsFixed(weirdo),2);
%         end
    end
end

%% Check if passive tuning was done
MyTuningTrials = []; TuningTrialSequence = []; PassiveReplayTraces = []; TuningParams = [];
[TuningFile] = WhereTuningFile(FilePaths,MyFileName);
if ~isempty(TuningFile)
    [MyTuningTrials, TuningTrialSequence, PassiveReplayTraces, TuningParams] = ParseTuningSession(TuningFile);
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
            GetOepsAuxChannels(mySortingdir, Trials.TimeStamps, MyTuningTrials, TuningTrialSequence, 'KiloSorted', 1); % send 'ADC', 1 to also get analog aux data
end

if isempty(TTLs)
    [myephysdir] = WhereOEPSFile(MyFileName,FilePaths); % returns empty if no recording file was found
    if ~isempty(myephysdir)
        if size(myephysdir,1) == 1
            [TTLs,ReplayTTLs,TuningTTLs,~] = ...
                GetOepsAuxChannels(myephysdir, Trials.TimeStamps, MyTuningTrials, TuningTrialSequence); % send 'ADC', 1 to also get analog aux data
        else
            while isempty(TTLs) && ~isempty(myephysdir)
                [TTLs,ReplayTTLs,TuningTTLs,~] = ...
                    GetOepsAuxChannels(myephysdir(1,:), Trials.TimeStamps, MyTuningTrials, TuningTrialSequence);
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
    % fix problematic odor starts
    if ~isempty(problematic_odor_starts)
        OEPS_OdorStarts = TTLs.Trial(problematic_odor_starts,4);
        if ~any(abs(TrialInfo.OdorStart(problematic_odor_starts,1)-OEPS_OdorStarts)>=0.01)
            TrialInfo.OdorStart(problematic_odor_starts,2) = TrialInfo.OdorStart(problematic_odor_starts,1);
        else
            disp('something funky with computing odorstart from Lever trace and/or InRewardZone');
            keyboard;
%             if ~any(TrialInfo.OdorStart(InitiationsFixed(weirdo),2)>-1)
%                 % Initiation hold was larger than a second - that's why couldn't
%                 % compute it accurately from trial traces
%                 TrialInfo.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo.OdorStart(InitiationsFixed(weirdo),2);
%             end
        end
    end
end

if ~isempty(TuningTTLs)
    Tuningextras.sessionsettings = TuningParams;
    Tuningextras.sequence = TuningTrialSequence;
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
    else
        % try the sorting directory instead
        if exist(mySortingdir)
            FileLocations.Spikes = mySortingdir;
            SingleUnits = GetSingleUnits(mySortingdir);
            [SingleUnits] = Spikes2Trials(SingleUnits, TTLs.Trial(1:size(TrialInfo.TrialID,2),:), TuningTTLs);
        end
    end
end

%% Photometry
FTrace = [];
if ~isempty(TTLs)
    
end

%% Sniffs
% [SniffTS] = ReadThermistorData(MyFilePath); % in behavior timestamps
if any(ismember(DataTags,'thermistor'))
    [SniffTS] = Thermistor2Sniffs(MyData(:,[1 find(ismember(DataTags,'thermistor')) find(ismember(DataTags,'Motor'))]));
else
    SniffTS = [];
end

%% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
if ~any(abs(TTLs.Trial(1:numel(TrialInfo.Odor),2) - (TrialInfo.SessionTimestamps(:,2) + TrialStart_Ephys - TrialStart_behavior))>0.04)
    % factor to convert all behavior timestamps to match Ephys
    TimestampAdjust.ClosedLoop = TrialStart_Ephys - TrialStart_behavior; % add to behavior TS to convert to OEPS
else
    disp('clock drift in ephys and behavior files');
    keyboard;
end

[SniffTS_passive] = ReadThermistorData(TuningFile); % in behavior timestamps
% for the passive case, convert sniff timestamps to OEPS base
% check first for clock drifts - compare trial starts between OEPS and
% matlab
if ~any(abs((TuningTTLs(:,1) - TuningTTLs(1,1)) - (TuningTTLs(:,9) - TuningTTLs(1,9)))>0.04)
    TimestampAdjust.Passive = TuningTTLs(1,1) - TuningTTLs(1,9);
else
    disp('trial start mismatch in ephys and behavior tuning files');
    keyboard;
end

%% Saving stuff in one place
if ~exist(fileparts(savepath),'dir')
    mkdir(fileparts(savepath));
end
save(savepath, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', 'FileLocations', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits', 'Tuningextras', 'SniffTS', 'SniffTS_passive', 'TimestampAdjust','-append');
    
end