
WhereSession = '/mnt/data/Sorted/Q5/2022-11-22_16-21-07';
%WhereSession = '/mnt/data/Sorted/Q9/2022-11-19_17-14-37';
WhereSession = '/mnt/data/Sorted/Q8/2022-12-09_14-37-11';

myKsDir = WhereSession;
if exist(fullfile(myKsDir,'kilosort4')) == 7
    KS4SingleUnits = GetSingleUnitsKS4(fullfile(myKsDir,'kilosort4'), 3);
end

[~,MouseName] = fileparts(fileparts(myKsDir));
[~,filename] = fileparts(myKsDir);
filename = [MouseName,'_',regexprep(filename(1,1:10),'-',''),'_r0_processed.mat'];
savepath = '/mnt/data/';
if exist('KS4SingleUnits','var')
    save(fullfile(savepath,'forWDW',filename),"KS4SingleUnits",'-append');
end