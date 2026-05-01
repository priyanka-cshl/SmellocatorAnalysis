Mice = {'Q3', 'Q4', 'Q5', 'Q8', 'Q9','S6', 'S7', 'S12', 'O3'};
for m = 1:5
    SortingMain = fullfile('/mnt/storage/Sorted',Mice{m});
%     %SortingMain = fullfile('/mnt/data/EarlySorted',Mice{m});
%     SortingMain = fullfile('/media/priyanka/ABC-ntfs/EphysSorted',Mice{m});
    AllFolders = dir(SortingMain);
    for x = 1:size(AllFolders,1)
        if ~strcmp(AllFolders(x).name,'.')&&~strcmp(AllFolders(x).name,'..')...
                &&~contains(AllFolders(x).name,'*')&&~contains(AllFolders(x).name,'~')
        thisDir = fullfile(SortingMain,AllFolders(x).name,'kilosort4');
        disp(thisDir);
        try
            QuickSniffTTLMapper_v2(thisDir,1);
        catch
            keyboard;
        end
        HDMain = fullfile('/media/priyanka/ABC-ntfs/EphysSorted',Mice{m},AllFolders(x).name);
        copyfile(fullfile(SortingMain,AllFolders(x).name,'quickprocesssniffs.mat'), ...
            fullfile(HDMain,'quickprocesssniffs.mat'));
        end
    end
end