%MyFilePath = 'Q3_20221011_r0.mat';

Paths = WhichComputer();

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

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
FileLocations.Behavior = MyFilePath;
[FilePaths, MyFileName] = fileparts(MyFilePath);

%% Parse into trials
%[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[Trials,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

for itr = 1:12
    LeverHists{itr} = [];
end

figure;
%% trial wise
for trial = 1:size(TrialInfo.TrialID,2)
    subplot(4,3,TrialInfo.TargetZoneType(trial));
    hold on
    plot(Traces.Lever{trial}(501:TrialInfo.TimeIndices(trial,2)));    
    
    LeverHists{TrialInfo.TargetZoneType(trial)} = vertcat(...
        LeverHists{TrialInfo.TargetZoneType(trial)}, ...
        Traces.Lever{trial}(501:TrialInfo.TimeIndices(trial,2)) );
end

%%
figure('Name',MyFileName);
for itr = 1:12
    subplot(12,1,itr);
    histogram(LeverHists{itr},(0:0.25:5));
end
set(gcf,'Position',[2283 1 240 995]);
saveas(gcf,fullfile('/home/priyanka/Desktop/random',[MyFileName,'.png']));
close all;
