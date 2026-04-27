%% function to parse behavioral data from the mouse lever task
% stripped down version to obtain respiration data and basic trial info

function [Traces, TrialInfo, MyData, DataTags] = PreprocessSmellocatorDataLite(MyFilePath)

Paths = WhichComputer();

%% globals
% global MyFileName;
global SampleRate;
SampleRate = 500; % Samples/second
global startoffset;
startoffset = 1; % in seconds
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

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
FileLocations.Behavior = MyFilePath;
[FilePaths, MyFileName] = fileparts(MyFilePath);
disp(MyFileName);

%% Parse into trials
[Trials,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

% sanity check - did some guess work in CorrectMatlabSampleDrops to compute odor start - check if it made sense
problematic_odor_starts = [];
if ~isempty(InitiationsFixed)
    if any(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01)
        weirdo = find(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01);
        problematic_odor_starts = InitiationsFixed(weirdo);
    end
end

if ~isempty(problematic_odor_starts)
    keyboard;
end
    
end