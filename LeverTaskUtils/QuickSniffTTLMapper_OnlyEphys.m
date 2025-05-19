function [AllSniffs] = QuickSniffTTLMapper_OnlyEphys(myKsDir)
% myKsDir = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Ephys/Q9/2022-11-19_17-14-37/';
% location where the sorted output is - contains cluster info and a TTL
% file generated during Kilosort

%% 1: Where is the corresponding behavior/tuning file
%   - need this to load trials to calculate timestamp difference
%   - and then process respiration with corrected timestamps

% where to look
[Paths] = WhichComputer();
if strcmp(myKsDir(end),filesep)
    myKsDir = myKsDir(1:end-1);
end

% have sniffs already been processed
% if exist(fullfile(myKsDir,'quickprocesssniffs.mat'))
%     load (fullfile(myKsDir,'quickprocesssniffs.mat'));
%     return;
% end

[~, RespirationData, AllSniffs] = GetSniffsFrom_myauxfile(myKsDir);

save(fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs', 'RespirationData');

end
