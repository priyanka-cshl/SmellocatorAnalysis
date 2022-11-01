function [myephysdir] = WhereOEPSTuningFile(BehaviorFile,BehaviorFolder)
foo = strsplit(BehaviorFile,'_');
mousename = char(foo(1));
date = char(foo(2));
datetoken = [date(1:4),'-',date(5:6),'-',date(7:8)];
    
[Paths] = WhichComputer();
if exist(fullfile(Paths.Mapping.EphysRaw,mousename),'dir')
    root = fullfile(Paths.Mapping.EphysRaw,mousename);
else
    disp('Please specify the location of raw Ephys Data');
    keyboard;
end
myfolders = dir ([root,'/',datetoken,'*']);

if size(myfolders,1) == 0
    myephysdir = [];
    disp('no oeps file found');
elseif size(myfolders,1) == 1 % only one recording session found
    myephysdir = myfolders.name;
    myephysdir = fullfile(root,myephysdir);
    if ~isempty(dir(fullfile(myephysdir, sprintf('Record Node*'))))
        tempfolder = (dir(fullfile(myephysdir, sprintf('Record Node*'))));
        myephysdir = fullfile(myephysdir,tempfolder.name);
    end
else
    % get sizes of all behavior tuning files
    AllFiles = dir(fullfile(BehaviorFolder,['*',date,'*o*.mat']));
    for i = 1:numel(AllFiles)
        fSize(i) = AllFiles(i).bytes;
    end
    nUnique = unique(fSize,'stable');
    thisFile = dir(fullfile(BehaviorFolder,[BehaviorFile,'.mat']));
    % same # of tuning and recording files
    if size(myfolders,1) == numel(nUnique)
        whichFolder = find(nUnique==thisFile.bytes);
        myephysdir = myfolders(whichFolder).name;
        myephysdir = fullfile(root,myephysdir);
        if ~isempty(dir(fullfile(myephysdir, sprintf('Record Node*'))))
            tempfolder = (dir(fullfile(myephysdir, sprintf('Record Node*'))));
            myephysdir = fullfile(myephysdir,tempfolder.name);
        end
    else
        disp('recording file vs. tuning file numbers are confusing');
        keyboard;
%         for i = 1:size(myfolders,1)
%             ephysdir = myfolders(i).name;
%             ephysdir = fullfile(root,ephysdir);
%             if ~isempty(dir(fullfile(ephysdir, sprintf('Record Node*'))))
%                 tempfolder = (dir(fullfile(ephysdir, sprintf('Record Node*'))));
%                 ephysdir = fullfile(ephysdir,tempfolder.name);
%             end
%             myephysdir(i,:) = ephysdir;
%         end
    end
end

end