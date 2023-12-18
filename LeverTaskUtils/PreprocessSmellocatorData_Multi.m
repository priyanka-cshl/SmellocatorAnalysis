%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = PreprocessSmellocatorData_Multi(MyFilePaths,overwriteflag)

if size(MyFilePaths,1) < 2
    disp('Input multiple closed loop sessions or use PreprocessSmellocatorData.m');
    return;
end

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
global NameTag
NameTag = '';

clear Traces TrialInfo Trials 

for sess = 1:size(MyFilePaths,1) % how many behavior files
    
    MyFilePath = MyFilePaths(sess,:);
    
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
    
    [MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
    FileLocations.Behavior(sess,:) = MyFilePath;
    [FilePaths, MyFileName] = fileparts(MyFilePath);
    disp(MyFileName);
    
    %% Parse into trials
    %[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
    [Trials_temp,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
    [MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
    [Traces_temp, TrialInfo_temp] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials_temp);
    
    % sanity check - did some guess work in CorrectMatlabSampleDrops to compute
    % odor start - check if it made sense
    if ~isempty(InitiationsFixed)
        if any(abs(diff(TrialInfo_temp.OdorStart(InitiationsFixed,:),1,2))>=0.01)
            weirdo = find(abs(diff(TrialInfo_temp.OdorStart(InitiationsFixed,:),1,2))>=0.01);
            if any(TrialInfo_temp.OdorStart(InitiationsFixed(weirdo),2)>-1)
                disp('something funky with computing odorstart from Lever trace');
                keyboard;
                TrialInfo_temp.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo_temp.OdorStart(InitiationsFixed(weirdo),2);
            else
                % Initiation hold was larger than a second - that's couldn't
                % compute it accurately from trial traces
                TrialInfo_temp.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo_temp.OdorStart(InitiationsFixed(weirdo),2);
            end
        end
    end
    
    %% append stuff across behavior sessions
    Traces{sess}    = Traces_temp;
    TrialInfo{sess} = TrialInfo_temp;
    Trials{sess}    = Trials_temp;
    NameTag         = [NameTag, '_', MyFileName(end-1:end)];
    
    %% get times for when files were created
    X = dir(MyFilePath);
    FileTS(sess,:)  = datevec(X.date);
    
    %% sniffs
    [SniffTS{sess}] = ReadThermistorData(MyFilePath); % in behavior timestamps

end

%% Check if passive tuning was done
MyTuningTrials = []; TuningTrialSequence = []; PassiveReplayTraces = []; TuningParams = [];
[TuningFile] = WhereTuningFile(FilePaths,MyFileName);

if size(TuningFile,1)>=1
    % more than one tuning file found
    [MyTuningTrials, TuningTrialSequence, PassiveReplayTraces, TuningParams, TuningFileTS] = ParseTuningSession_Multi(TuningFile);
    disp(['Found ', num2str(size(TuningFile,1)), ' Tuning Files']);
    % FileLocations.Tuning = TuningFile;
else
    disp('Error: No tuning file found');
    return;
end

%% get the order of closed-loop and passive sessions - r's vs. o's
global FileOrder
[~, order] = sortrows([FileTS; TuningFileTS]);
FileOrder = [(1:size(FileTS,1))'; -1*(1:size(TuningFileTS,1))'];
FileOrder = FileOrder(order);

%% Get info from the OEPS files if available
TTLs = []; ReplayTTLs = []; TuningTTLs = []; myephysdir = [];

% first check if event data has been processed during Sorting - this
% should be the default for any mid-session replays
[mySortingdir] = WhereOEPSTTLs(MyFileName,FilePaths); % returns empty if no recording file was found
if ~isempty(mySortingdir)
    % process the TTLs
    [TTLs,ReplayTTLs,BehaviorTrials,TuningTTLs,~] = ...
            GetOepsAuxChannels_Multi(mySortingdir, Trials, MyTuningTrials, TuningTrialSequence, 'KiloSorted', 1); % send 'ADC', 1 to also get analog aux data
end

if isempty(TTLs)
    disp('no matching recording file found');
else
    FileLocations.KiloSorted = mySortingdir;
    FileLocations.OEPS = myephysdir;
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
        [SingleUnits, TrialStretches] = Spikes2Trials_Multi(SingleUnits, TTLs.Trial, BehaviorTrials, TuningTTLs);   
    else
        % try the sorting directory instead
        if exist(mySortingdir)
            FileLocations.Spikes = mySortingdir;
            SingleUnits = GetSingleUnits(mySortingdir);
            [SingleUnits, TrialStretches] = Spikes2Trials_Multi(SingleUnits, TTLs.Trial, BehaviorTrials, TuningTTLs);   
        end
    end
end

%% Sniffs
for n = 1:size(TuningFile,1)
    if ~isempty(TuningTTLs{n})
        [SniffTS_passive{n}] = ReadThermistorData(TuningFile(n,:)); % in behavior timestamps
        % for the passive case, convert sniff timestamps to OEPS base
        % check first for clock drifts - compare trial starts between OEPS and
        % matlab
        
        if ~any(abs((TuningTTLs{n}(:,1) - TuningTTLs{n}(1,1)) - (TuningTTLs{n}(:,9) - TuningTTLs{n}(1,9)))>0.04)
            Passive_Timestamp_adjust(n) = TuningTTLs{n}(1,1) - TuningTTLs{n}(1,9);
        else
            disp('trial start mismatch in ephys and behavior tuning files');
            keyboard;
        end
    else
        Passive_Timestamp_adjust(n) = NaN;
    end
end

%% concatenate the multiple closed-loop sessions for easier processing later
[Traces, TrialInfo, SniffTS] = StitchBehavior(Traces,TrialInfo,SniffTS,BehaviorTrials);

%% Saving stuff in one place
savepath = fullfile(Paths.Grid.Behavior_processed,AnimalName,[regexprep(MyFileName,'_r[0-9]',NameTag),'_processed.mat']);

if ~exist(fileparts(savepath),'dir')
    mkdir(fileparts(savepath));
end
save(savepath, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', 'FileLocations', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits', 'Tuningextras', 'SniffTS', 'SniffTS_passive', 'Passive_Timestamp_adjust', ...
               'TrialStretches');
    
end