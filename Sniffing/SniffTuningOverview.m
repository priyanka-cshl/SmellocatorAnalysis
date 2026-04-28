%% script for loading the quick sorted sessions for early Smellocator data
%  and ephys

%% Path handling
SessionName = 'Q3_20220906';
if ~contains(SessionName,'.mat')
    if contains(mfilename('fullpath'),'opt') % Linux desktops
        basePath = '/opt/';
        FileSplit = strsplit(SessionName,'_');
        if regexp(FileSplit{1},'Q\d') 
            MouseName = FileSplit{1};
            if regexp(FileSplit{2},'\d*') % date
                DateStr = FileSplit{2};
            end
        end

        MouseName = SessionName(1:strfind(SessionName,'_')-1);
        filesearch = dir(fullfile('/mnt/storage/Sorted',MouseName,[Ses])
        sessionPath = fullfile('/mnt/storage/Sorted', SessionName);
    end
end


% SessionName = 'APCB_20220311_r0_processed.mat';
% SessionName = 'Q3_20221026_r0_processed.mat';
% SessionName = 'Q4_20221112_r0_processed.mat';
%SessionName = 'Q5_20221122_r0_processed.mat';
%SessionName = 'Q8_20221209_r0_processed.mat';
% SessionName = 'Q9_20221119_r0_processed.mat';
SessionName = 'T3_20250508_r0_processed.mat';
%SessionName = 'T2_20250508_r0_processed.mat';

%% Housekeeping
% Some path handling
if exist('/home/wolf/Documents/')
    basePath = '/home/wolf/Documents/';    
    sessionPath = fullfile('/home/wolf/Documents/SmellocatorGLM/data/sessions',...
                  strrep(SessionName,'.mat',''), ...
                  'data.mat');
elseif contains(mfilename('fullpath'),'opt') % Linux desktops
    basePath = '/opt/';
    sessionPath = fullfile('/mnt/data/forWDW', SessionName);
else % macbook
    basePath = '/Users/Priyanka/Desktop/github_local';
    sessionPath = fullfile(...
                  '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/forWDW', ...
                  SessionName);
end
addpath(fullfile(basePath,'MatlabUtils'));

%%  Loading the processed data 
data = load(sessionPath);

if strfind(SessionName,'Q')
    [AllSniffs, SniffColumnInfo] = GetAllSniffs(sessionPath,'whichsensor',1);
    SniffStarts = [];
    SniffStarts(:,3) = AllSniffs(1:end-2,1);
    SniffStarts(:,4) = AllSniffs(2:end-1,1);
    SniffStarts(:,5) = AllSniffs(3:end,1);
    SniffStarts(:,6) = (SniffStarts(:,5) - SniffStarts(:,4)).^-1; % this sniff freq
    SniffStarts(:,7) = (SniffStarts(:,4) - SniffStarts(:,3)).^-1; % prev sniff freq
    SniffStarts(:,8) = AllSniffs(2:end-1,4); % manifold state
    SniffStarts(:,9) = AllSniffs(2:end-1,5); % odor state
    SniffStarts(:,10) = AllSniffs(2:end-1,8); % session phase
    SniffStarts(:,11) = AllSniffs(2:end-1,6); % odor location

    SniffStarts(SniffStarts(:,10)~=1,:) = [];
    SniffStarts(SniffStarts(:,8)>0 & SniffStarts(:,8)<1,:) = [];
    SniffStarts(SniffStarts(:,9)<0,:) = [];

    % only keep the ITI sniffs
    f = find(SniffStarts(:,8)~=0 | SniffStarts(:,9)~=0);
    SniffStarts(f,:) = [];
else

    %% get sniff times and sniff freq
    timeResolution = mode(diff(data.TracesOut.Timestamps{1}));
    SniffStarts = find(diff(data.TracesOut.SniffsDigitized{1})==1);
    SniffStarts(:,3) = circshift(SniffStarts,-1);
    SniffStarts(:,2) = SniffStarts(:,1);
    SniffStarts(:,1) = circshift(SniffStarts(:,2),1);
    SniffStarts(1,1) = nan;
    SniffStarts(end,3) = nan;
    SniffStarts([1 end],:) = [];
    SniffStarts(:,4) = data.TracesOut.Timestamps{1}(SniffStarts(:,2)); % timestamp
    SniffStarts(:,5) = data.TracesOut.Timestamps{1}(SniffStarts(:,3)); % timestamp
    SniffStarts(:,6) = ((SniffStarts(:,3) - SniffStarts(:,2))*timeResolution).^-1; % this sniff freq
    SniffStarts(:,7) = ((SniffStarts(:,2) - SniffStarts(:,1))*timeResolution).^-1; % prev sniff freq
    SniffStarts(:,8:10) = 0;

end

maxsniffduration = 0.7;
f = find(SniffStarts(:,5)-SniffStarts(:,4)>maxsniffduration);
SniffStarts(f,5) = SniffStarts(f,4) + maxsniffduration;

%% plot units from KS4 or KS1
%get clusterIDs and tetrode entries for all units
if isfield(data,'KS4SingleUnits')
    myUnits = [[data.KS4SingleUnits.id]' [data.KS4SingleUnits.tetrode]'];
    unittag = 'KS4SingleUnits';
elseif isfield(data,'SingleUnits')
    myUnits = [[data.SingleUnits.id]' [data.SingleUnits.tetrode]'];
    unittag = 'SingleUnits';
else
    keyboard;
end

myUnits(:,3) = floor(myUnits(:,2));
myUnits(:,3) = floor((1:size(myUnits,1))/10)+1;
% myUnits(:,3) = 1;

%%
figure;
for snifforder = 3 %1:3
    switch snifforder
        case 1
            SniffStarts = sortrows(SniffStarts, [10 8 9 1], 'descend');
        case 2
            % we can resort SniffStarts by current sniff frequency
            SniffStarts = sortrows(SniffStarts, [10 8 9 6], 'descend');
        case 3
            % we can resort SniffStarts by current sniff frequency
            SniffStarts = sortrows(SniffStarts, [10 8 9 7], 'descend');
    end

    for whichtetrode = 1:16 %16
        f = find(myUnits(:,3)==whichtetrode);
%         figure; %(whichtetrode); 
        %subplot(5,4,whichtetrode);
        subplot(5,2,whichtetrode);
        hold on
%        subplot(1,3,whichplot); hold on
%         whichplot = snifforder;
%         subplot(1,3,whichplot); hold on
        %whichplot = (s-1)*3 + snifforder;
        %subplot(3,3,whichplot); hold on
        for i = 1:numel(f)
            whichunit = f(i);
            spikelist = thisUnitSniffPlot(data.(unittag)(whichunit).spikes,SniffStarts(:,4:5));
            offset = (i-1)*maxsniffduration;
            %plot(spikelist(:,2)+offset,spikelist(:,1),'.');
            plot(spikelist(:,2)+offset,spikelist(:,3),'.','MarkerSize',2);
        end
        set(gca, 'XLim',[0 11])
    end
    
end

function [spikelist] = thisUnitSniffPlot(spikes,snifftimes)
spikelist = [];
for s = 1:size(snifftimes,1)-1
    t = [snifftimes(s,1)  snifftimes(s,2)];
    thisSniffSpikes = spikes(spikes>=t(1)&spikes<t(2));
    thisSniffSpikes(:,2) = thisSniffSpikes - snifftimes(s,1);
    thisSniffSpikes(:,3) = s + 0*thisSniffSpikes(:,1);
    spikelist = vertcat(spikelist, thisSniffSpikes);
end
end





