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
if ~exist(fullfile(myKsDir,'quickprocesssniffs.mat'))
    [~, RespirationData, AllSniffs] = GetSniffsFrom_myauxfile(myKsDir);
else
    load(fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs', 'RespirationData');
end

% add odor TTL info if available
if exist(fullfile(myKsDir,'quickprocessOdorTTLs.mat'))
    load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings');
    
    % add the column infos to the sniffs
    for i = 1:size(TTLs.Trial,1) % everty trial
        t = TTLs.Trial(i,7:10); % odor1 start, stop, purge start, stop
        t(end+1) = t(end) + TTLs.Trial(i,9)-TTLs.Trial(i,8); % add a post-purge period the same as the first pulse
        TS = [t(1:end-1)' t(2:end)' [1 3 2 4]'];
        for j = 1:size(TS,1)
            whichsniffs = find( (AllSniffs(:,1)>=TS(j,1)) & (AllSniffs(:,1)<TS(j,2)) );
            AllSniffs(whichsniffs,5) = TTLs.Trial(i,4); % odor identity
            AllSniffs(whichsniffs,6) = TS(j,3);
        end
    end
    
end

save(fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs', 'RespirationData');

end
