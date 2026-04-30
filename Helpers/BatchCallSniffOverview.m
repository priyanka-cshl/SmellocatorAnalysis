Mice = {'Q3', 'Q4', 'Q5', 'Q8', 'Q9'};
for m = 1:5
    %figure;
    %SortingMain = fullfile('/mnt/storage/Sorted',Mice{m});
    %SortingMain = fullfile('/mnt/data/EarlySorted',Mice{m});
    SortingMain = fullfile('/media/priyanka/ABC-ntfs/EphysSorted',Mice{m});
    AllFolders = dir(SortingMain);
    y = 0;
    for x = 1:size(AllFolders,1)
        if ~strcmp(AllFolders(x).name,'.')&&~strcmp(AllFolders(x).name,'..')...
                &&~contains(AllFolders(x).name,'*')&&~contains(AllFolders(x).name,'~')
        thisDir = fullfile(SortingMain,AllFolders(x).name,'kilosort4');
        disp(thisDir);
        y = y + 1;
        %% Sniff distributions
        sessionPath = fullfile(SortingMain,AllFolders(x).name);
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
        end

        %% Keep only a subset of sniffs
        %subplot(1,5,y); hold on
        for whichsniffs = -2:1:3
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
            H(:,whichsniffs+3,y) = histcounts(ValidSniffs(:,6),[0:1:12]);
            %histogram(ValidSniffs(:,6));
        end
        %plot(H);

        %%
        %SniffTuningOverview(thisDir);
        

%         try
%             QuickSniffTTLMapper_v2(thisDir,1);
%         catch
%             keyboard;
%         end
%         HDMain = fullfile('/media/priyanka/ABC-ntfs/EphysSorted',Mice{m},AllFolders(x).name);
%         copyfile(fullfile(SortingMain,AllFolders(x).name,'quickprocesssniffs.mat'), ...
%             fullfile(HDMain,'quickprocesssniffs.mat'));
        end

    end
    figure;
    for c = 1:y
        subplot(5,1,c); hold on
        plot(H(:,:,c));
    end
end


% Q3 - 2nd
% Q4 - 1st, last, 2nd also good
% Q5 - 1st, last
% Q8 - 1st
% Q9 - 1st