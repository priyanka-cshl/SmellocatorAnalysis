
%% Script for running behavior related analysis on multiple sessions in a loop
function [] = LoopThroughSessions(MouseName)

%% File selection
[Paths] = WhichComputer(); % load rig specific paths etc
if nargin
    DataRoot = fullfile(Paths.Grid.Behavior,MouseName);
else
    DataRoot = uigetdir([],'Select the folder containing behavior sessions');
    [~,MouseName] = fileparts(DataRoot); %#ok<ASGLU>
    %[FileNames,FilePaths] = uigetfile('*_r*.mat','choose one or more session files','MultiSelect','on',DataRoot);
end

AllFiles = dir(fullfile(DataRoot,'*_r*.mat'));

for i = 1:size(AllFiles.name,2) % For each file
    PreprocessSmellocatorData(fullfile(DataRoot,AllFiles(i).name));
end
