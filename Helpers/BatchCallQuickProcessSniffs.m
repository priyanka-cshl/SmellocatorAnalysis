Mice = {'Q3', 'Q4', 'Q5', 'Q8', 'Q9','S6', 'S7', 'S12', 'O3'};
Mice = {'E2'};
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
% for batch copying of quick process sniffs and odor TTLs to Hard drives
% for at home processing of sniff times or analysis
% quickprocesssniffs.mat can have SingleUnit info already

Mice = {'E2','E3','E6'};
BasePath    = '/mnt/data/Sorted'; % on Catalina
HDBasePath  = '/media/priyanka/ABC-ntfs/EphysSorted';
for m = 1:numel(Mice)
    SortingMain = fullfile(BasePath,Mice{m});
    AllFolders = dir(SortingMain);
    for x = 1:size(AllFolders,1)
        if ~strcmp(AllFolders(x).name,'.')&&~strcmp(AllFolders(x).name,'..')
            thisDir = fullfile(SortingMain,AllFolders(x).name);
            disp(thisDir);
            HDpath = fullfile(HDBasePath,Mice{m},AllFolders(x).name);
            if ~exist(HDpath, 'dir')
                mkdir(HDpath)
            end
            copyfile(fullfile(thisDir,'quickprocesssniffs.mat'), ...
                    fullfile(HDpath,'quickprocesssniffs.mat'));
            copyfile(fullfile(thisDir,'quickprocessOdorTTLs.mat'), ...
                    fullfile(HDpath,'quickprocessOdorTTLs.mat'));
        end
    end
end