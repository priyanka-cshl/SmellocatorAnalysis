%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = PreprocessCIDData(MyFilePath,overwriteflag)

if nargin<2
    overwriteflag = 0;
end

%% Add relevant repositories
Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));
addpath(genpath(fullfile(Paths.Code,'open-ephys-matlab-tools'))); % for the new OEPS GUI: https://github.com/open-ephys/open-ephys-matlab-tools (use commit 10-04-22)
addpath(genpath(fullfile(Paths.Code,'MatlabUtils')));
addpath(genpath(fullfile(Paths.Code,'npy-matlab/'))); % path to npy-matlab scripts

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
savepath = fullfile(Paths.CID.Processed,AnimalName,[MyFileName,'_cid-processed.mat']);
if exist(savepath)
    if ~overwriteflag 
        reply = input('This session has already been processed. \nDo you want to overwrite? Y/N [Y]: ','s');
        if strcmp(reply,'N')
            return;
        end
    end
    load(savepath,'SniffData','TTLs');
end

%% Have sniff data and TTLs been extracted already? 
if ~exist('SniffData','var')
    %rawdata = fullfile(Paths.CID.RawFiles,AnimalName,'raw',MyFileName);
    rawdata = fullfile(FilePaths,'raw',MyFileName);
    GetOepsTTLsAndRespiration_CID(rawdata);
    load(savepath,'SniffData','TTLs');
end

%% trial and odor traces
Timestamps = SniffData(:,1);
StimTrace  = zeros(numel(Timestamps),4); % [trialOn odorIdentity odorIntensity]
for i = 1:size(TTLs.Trial,1) % every trial
    % trial on
    [~,idx1] = min(abs(Timestamps-TTLs.Trial(i,1))); % trial start
    [~,idx2] = min(abs(Timestamps-TTLs.Trial(i,2))); % trial stop
    StimTrace(idx1:idx2,1) = 1;
    % odor on
    [~,idx1] = min(abs(Timestamps-TTLs.Trial(i,7))); % stimulus start
    [~,idx2] = min(abs(Timestamps-TTLs.Trial(i,8))); % stimulus stop
    if TTLs.Trial(i,4)
        StimTrace(idx1:idx2,2) = TTLs.Trial(i,4); % stimulus identity
    else
        StimTrace(idx1:idx2,2) = max(TTLs.Trial(:,4)) + 1;
    end
    StimTrace(idx1:idx2,3) = TTLs.Trial(i,5); % stimulus intensity
end
TracesOut.Timestamps{1}     = Timestamps;
TracesOut.Trial{1}  = StimTrace(:,1);
TracesOut.Odor{1}   = StimTrace(:,2);
TracesOut.OdorIntensity{1}   = StimTrace(:,3);

%% Sniffs
TracesOut.Sniffs{1}         = SniffData(:,2);
TracesOut.SniffsFiltered{1} = SniffData(:,3);

% add a digital sniff trace
load(savepath,'CuratedSniffTimestamps');
if exist('CuratedSniffTimestamps','var')
    if size(CuratedSniffTimestamps,2) < 10
        CuratedSniffTimestamps(:,10) = 0;
    end
    DigitalSniffs = TracesOut.SniffsFiltered{1}*0;
    for n = 1:size(CuratedSniffTimestamps)
        idx = CuratedSniffTimestamps(n,8:9);
        if CuratedSniffTimestamps(n,10) == -1
            DigitalSniffs(idx(1):idx(2)) = -1;
        else
            DigitalSniffs(idx(1):idx(2)) = 1;
        end
    end
    TracesOut.SniffsDigitized{1} = DigitalSniffs;
end

%% Get spikes - label spikes by trials
SingleUnits = GetSingleUnits(MyFilePath);

%% append to the existing processed file
save(savepath, 'TracesOut','SingleUnits','-append');
    
end