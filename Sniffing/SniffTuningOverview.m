%% script for loading the quick sorted sessions for early Smellocator data
%  and ephys

%% Path handling
% SessionName = 'Q3_20220906';
SessionName = 'Q5_20220907';
% SessionName = 'Q8_20220928';
% SessionName = 'Q4_20220906';
SessionName = 'Q9_20220927';
SessionName = 'Q9_20221001';
if ~contains(SessionName,'.mat')
    if contains(mfilename('fullpath'),'opt') % Linux desktops
        basePath = '/opt/';
        FileSplit = strsplit(SessionName,'_');
        if regexp(FileSplit{1},'Q\d') 
            MouseName = FileSplit{1};
            if regexp(FileSplit{2},'\d*') % date
                DateStr = FileSplit{2};
                DateToken = [DateStr(1:4),'-',DateStr(5:6),'-',DateStr(7:8)];

                % find the ephys dir
                EphysDir = dir(fullfile('/mnt/storage/Sorted',MouseName,[DateToken,'*']));
                sessionPath = fullfile(EphysDir.folder,EphysDir.name);
            end
        end
    end
end

%% Housekeeping
addpath(fullfile(basePath,'MatlabUtils'));

%% Load the sorted units from KS4
if exist(fullfile(sessionPath, 'kilosort4'),'dir')
    AllUnits = GetSingleUnits(fullfile(sessionPath, 'kilosort4'));
    unittag = 'KS4SingleUnits';
else
    AllUnits = GetSingleUnits(sessionPath);
    unittag = 'SingleUnits';
end
myUnits = [[AllUnits.id]' [AllUnits.tetrode]' [AllUnits.quality]'];
myUnits(:,4) = 1:size(myUnits,1);

if strcmp(unittag,'KS4SingleUnits') && isempty(dir(fullfile(sessionPath,'kilosort4','cluster_info*')))
    % session wasn't curated in phy, keep only 'good' units
    myUnits(find(myUnits(:,3)~=2),:) = [];
end
% myUnits(:,3) = floor(myUnits(:,2));
% myUnits(:,3) = floor((1:size(myUnits,1))/10)+1;

%% Load the sniffs
if exist(fullfile(sessionPath, 'quickprocesssniffs.mat'))
    RespData = load(fullfile(sessionPath, 'quickprocesssniffs.mat'));
    % pass info into SniffStarts
    SniffStarts = [];
    SniffStarts(:,3) = RespData.AllSniffs(1:end-2,1); % TS of prev inh
    SniffStarts(:,4) = RespData.AllSniffs(2:end-1,1); % TS of this inh
    SniffStarts(:,5) = RespData.AllSniffs(3:end,1);   % TS of next inhalation

    SniffStarts(:,1) = RespData.AllSniffs(2:end-1,1); % TS: start of this inh
    SniffStarts(:,2) = RespData.AllSniffs(2:end-1,2); % TS: end of this inh

    SniffStarts(:,6) = (SniffStarts(:,5) - SniffStarts(:,4)).^-1; % this sniff freq
    SniffStarts(:,7) = (SniffStarts(:,4) - SniffStarts(:,3)).^-1; % prev sniff freq

    SniffStarts(:,11) = RespData.AllSniffs(2:end-1,4); % odor location
    SniffStarts(:,10) = 1; % session phase

    % get manifold and odor state from the ephys TTLs
    if isfield(RespData,'Traces')
        TS = RespData.Traces.Timestamps{1};
        Manifold = RespData.Traces.Manifold{1};
        Odor = RespData.Traces.Odor{1};
        for i = 1:size(SniffStarts,1)
            x1 = find(TS>=SniffStarts(i,1),1,"first");
            x2 = find(TS>=SniffStarts(i,2),1,"first");
            SniffStarts(i,8) = mode(Manifold(x1:x2)); % manifold state
            SniffStarts(i,9) = mode(Odor(x1:x2)); % manifold state
        end
    end
    
    % SniffStarts(SniffStarts(:,9)<0,:) = [];
else
    keyboard;
end

%% Keep only a subset of sniffs
whichsniffs = -1;
switch whichsniffs
    case -1 % only keep the ITI sniffs : manifold off, air on
        f = (find(SniffStarts(:,8)==0 & SniffStarts(:,9)==0));
    case -2 % only keep the no air sniffs : manifold off, air off
        f = (find(SniffStarts(:,8)==0 & SniffStarts(:,9)==-1));
    case 0 % manifold ON, air ON
        f = (find(SniffStarts(:,8)==1 & SniffStarts(:,9)==0));
    case 1 % manifold ON, Odor1 ON
        f = (find(SniffStarts(:,8)==1 & SniffStarts(:,9)==1));
    case 2 % manifold ON, Odor2 ON
        f = (find(SniffStarts(:,8)==1 & SniffStarts(:,9)==2));
    case 3 % manifold ON, Odor3 ON
        f = (find(SniffStarts(:,8)==1 & SniffStarts(:,9)==3));
    case 4 % Manifold ON
        f = (find(SniffStarts(:,8)==1 & SniffStarts(:,9)>=0));
end

ValidSniffs = SniffStarts(f,:);

%% Plot sniff tuned rasters
nTetrodes = max(myUnits(:,2));
nUnits = size(myUnits,1);
[~,foo,~] = unique(floor(myUnits(:,2)));
unitsPerTetrode = circshift(foo(:,1),-1) - foo(:,1);
unitsPerTetrode(end) = nUnits - sum(unitsPerTetrode(1:end-1));
maxsniffduration = 0.7;

figure;
for snifforder = 3 %1:3
    switch snifforder
        case 1 % manifold state, occurence
            ValidSniffs = sortrows(ValidSniffs, [10 8 9 1], 'descend');
        case 2
            % we can resort ValidSniffs by current sniff frequency
            ValidSniffs = sortrows(ValidSniffs, [10 8 9 6], 'descend');
        case 3
            % we can resort ValidSniffs by previous sniff frequency
            ValidSniffs = sortrows(ValidSniffs, [10 8 9 7], 'descend');
    end

    for whichtetrode = 1:nTetrodes
        f = find(floor(myUnits(:,2))==whichtetrode);
        subplot(5,2,whichtetrode);
        hold on
        for i = 1:numel(f)
            whichunit = myUnits(f(i),4);
            spikelist = thisUnitSniffPlot(AllUnits(whichunit).spikes,ValidSniffs(:,4:5));
            % ignore any spikes beyond the max sniffduration
            spikelist(spikelist(:,2)>maxsniffduration,:) = [];
            offset = (i-1)*maxsniffduration;
            %plot(spikelist(:,2)+offset,spikelist(:,1),'.');


            plot(offset+(ValidSniffs(:,4)*0),1:size(ValidSniffs,1),'k');
            if snifforder == 2
                plot(offset+(ValidSniffs(:,5)-ValidSniffs(:,4)),1:size(ValidSniffs,1),'k');
            end

            plot(spikelist(:,2)+offset,spikelist(:,3),'.','MarkerSize',2);
        end
        set(gca, 'XLim',[0 max(unitsPerTetrode)*maxsniffduration]);
    end
    
end

%%
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





