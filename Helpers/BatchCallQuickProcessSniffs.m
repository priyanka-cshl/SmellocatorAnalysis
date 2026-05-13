Mice = {'Q3', 'Q4', 'Q5', 'Q8', 'Q9','S6', 'S7', 'S12', 'O3'};

%%
Mice = {'APC2'};
copytoHD = 0;
for m = 1:numel(Mice)
    SortingMain = fullfile('/mnt/data/Sorted',Mice{m});
    %   SortingMain = fullfile('/mnt/storage/Sorted',Mice{m});
    %     %SortingMain = fullfile('/mnt/data/EarlySorted',Mice{m});
    %     SortingMain = fullfile('/media/priyanka/ABC-ntfs/EphysSorted',Mice{m});
    AllFolders = dir(SortingMain);
    for x = 1:size(AllFolders,1)
        if ~strcmp(AllFolders(x).name,'.')&&~strcmp(AllFolders(x).name,'..')...
                &&~contains(AllFolders(x).name,'*')&&~contains(AllFolders(x).name,'~')
            thisDir = fullfile(SortingMain,AllFolders(x).name,'kilosort4');
            disp(thisDir);
            try
                %QuickSniffTTLMapper_v2(thisDir,1);
                myKsDir = fullfile(SortingMain,AllFolders(x).name);
                makeTTLs_CID(myKsDir); 
                QuickSniffTTLMapper_OnlyEphys(myKsDir,1);
            catch
                keyboard;
            end
            if copytoHD
                HDMain = fullfile('/media/priyanka/ABC-ntfs/EphysSorted',Mice{m},AllFolders(x).name);
                copyfile(fullfile(SortingMain,AllFolders(x).name,'quickprocesssniffs.mat'), ...
                    fullfile(HDMain,'quickprocesssniffs.mat'));
            end
        end
    end
end


%% 
% for batch processing of odor TTLs and sniff data for CID sessions

Mice = {'D4'};
BasePath    = '/mnt/data/Sorted'; % on Catalina
HDBasePath  = '/media/priyanka/ABC-ntfs/EphysSorted';
copytoHD    = 0; 
for m = 1:numel(Mice)
    SortingMain = fullfile(BasePath,Mice{m});
    AllFolders = dir(SortingMain);
    for x = 6:size(AllFolders,1)
        if ~strcmp(AllFolders(x).name,'.')&&~strcmp(AllFolders(x).name,'..')
            thisDir = fullfile(SortingMain,AllFolders(x).name);
            disp(thisDir);

            % get sniff data and odor TTLs
            makeTTLs_CID(thisDir); 
            QuickSniffTTLMapper_OnlyEphys(thisDir,1);
            
            if copytoHD
                HDpath = fullfile(HDBasePath,Mice{m},AllFolders(x).name);
                copyfile(fullfile(thisDir,'quickprocesssniffs.mat'), ...
                    fullfile(HDMain,'quickprocesssniffs.mat'));
                copyfile(fullfile(thisDir,'quickprocessOdorTTLs.mat'), ...
                    fullfile(HDMain,'quickprocessOdorTTLs.mat'));
            end
        end
    end
end

%% 
% for batch copying of quick process sniffs and odor TTLs to Hard drives
% for at home processing of sniff times or analysis
% quickprocesssniffs.mat can have SingleUnit info already

Mice = {'D4'};
BasePath    = '/mnt/data/Sorted'; % on Catalina
HDBasePath  = '/media/priyanka/ABC-ntfs/EphysSorted';
for m = 1:numel(Mice)
    SortingMain = fullfile(BasePath,Mice{m});
    AllFolders = dir(SortingMain);
    for x = 6:size(AllFolders,1)
        if ~strcmp(AllFolders(x).name,'.')&&~strcmp(AllFolders(x).name,'..')
            thisDir = fullfile(SortingMain,AllFolders(x).name);
            disp(thisDir);
            HDMain = fullfile(HDBasePath,Mice{m},AllFolders(x).name);
            if ~exist(HDMain,'dir')
                mkdir(HDMain);
            end
            copyfile(fullfile(thisDir,'quickprocesssniffs.mat'), ...
                    fullfile(HDMain,'quickprocesssniffs.mat'));
            copyfile(fullfile(thisDir,'quickprocessOdorTTLs.mat'), ...
                    fullfile(HDMain,'quickprocessOdorTTLs.mat'));
        end
    end
end

%% get unit summary and protocol summary of sessions
Mice = {'APC1'};
BasePath    = '/mnt/data/Sorted'; % on Catalina 
% BasePath  = '/media/priyanka/ABC-ntfs/EphysSorted';
for m = 1:numel(Mice)
    SortingMain = fullfile(BasePath,Mice{m});
    AllFolders = dir(SortingMain);
    for x = 1:size(AllFolders,1)
        if ~strcmp(AllFolders(x).name,'.')&&~strcmp(AllFolders(x).name,'..')
            myKsDir = fullfile(SortingMain,AllFolders(x).name);
            disp(myKsDir);
            load(fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs','KS4Units'); % sniff times

            if exist('KS4Units') && isempty(dir(fullfile(myKsDir,'kilosort4','cluster_info*')))
                myUnits = [[KS4Units.id]' [KS4Units.tetrode]' [KS4Units.quality]'];
                totalunits = size(myUnits,1);
                goodUnits = numel(find(myUnits(:,3)==2));
                disp(['found ',num2str(totalunits),' units; ', num2str(goodUnits),' good units']);
            end

            load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings'); % odorTTLs
            disp(StimSettings.SessionType);
            disp(['found ',num2str(size(TTLs.Trial,1)),' trials; ', num2str(numel(StimSettings.Odors)*StimSettings.Reps),' trials in stim file']);
        end
    end
end
