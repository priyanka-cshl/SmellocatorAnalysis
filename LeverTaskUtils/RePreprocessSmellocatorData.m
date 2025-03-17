%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = RePreprocessSmellocatorData(MyFilePath,newname)

%% Add relevant repositories
Paths = WhichComputer();
% addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));
% addpath(genpath(fullfile(Paths.Code,'open-ephys-matlab-tools'))); % for the new OEPS GUI: https://github.com/open-ephys/open-ephys-matlab-tools (use commit 10-04-22)
% addpath(genpath(fullfile(Paths.Code,'MatlabUtils')));
% addpath(genpath(fullfile(Paths.Code,'npy-matlab/'))); % path to npy-matlab scripts

foo = regexp(MyFilePath,'_','split');
AnimalName = foo{1};
MyFilePath = fullfile(Paths.Grid.Behavior,AnimalName,MyFilePath);
[FilePaths, MyFileName] = fileparts(MyFilePath); %#ok<*ASGLU>

%% check if the preprocessed version already exists - locally or on the server
oldpath = fullfile(Paths.Grid.Behavior_processed,AnimalName,[MyFileName,'_processed.mat']);
load(oldpath);

savepath = regexprep(oldpath,AnimalName,newname);
if ~exist(savepath)
    copyfile(oldpath,savepath);
end
%% Get spikes - label spikes by trials
SingleUnits = [];
% get directory to choose sorting units
[~,sortingdir] = fileparts(FileLocations.Spikes);
newsortingdir = uigetdir(fullfile(Paths.Local.Ephys_processed,AnimalName,sortingdir), 'Pick new sorting location');

if newsortingdir
    myspikesdir = newsortingdir;
    FileLocations.KiloSorted = myspikesdir;
    FileLocations.Spikes = myspikesdir;
    SingleUnits = GetSingleUnits(myspikesdir);
    [SingleUnits] = Spikes2Trials(SingleUnits, TTLs.Trial(1:size(TrialInfo.TrialID,2),:), TuningTTLs);
end


%% Saving stuff in one place
save(savepath, 'FileLocations', 'SingleUnits','-append');

end