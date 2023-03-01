function [mySortingdir] = WhereOEPSTTLs(BehaviorFile,BehaviorFolder)
foo = strsplit(BehaviorFile,'_');
mousename = char(foo(1));
date = char(foo(2));
datetoken = [date(1:4),'-',date(5:6),'-',date(7:8)];
    
[Paths] = WhichComputer();

% check if the Sorting has already been done, and if yes, TTLs have been
% pulled out

if exist(fullfile(Paths.Grid.Ephys_processed,mousename),'dir')
    root = fullfile(Paths.Grid.Ephys_processed,mousename);
elseif exist(fullfile(Paths.Local.Ephys_processed,mousename),'dir')
    root = fullfile(Paths.Local.Ephys_processed,mousename);
end

myfolders = dir ([root,'/',datetoken,'*']);

if size(myfolders,1) == 0
    mySortingdir = [];
    disp('no oeps file found');
elseif size(myfolders,1) == 1 % only one recording session found
    mySortingdir = myfolders.name;
    mySortingdir = fullfile(root,mySortingdir);
else
    % same # of behavior and recording files
%     if numel(dir(fullfile(BehaviorFolder,['*',date,'*r*.mat']))) == size(myfolders,1)
%         sessionname = char(foo(3));
%         session_num = str2num(sessionname(2));
%         mySortingdir = myfolders(session_num+1).name;
%         mySortingdir = fullfile(root,mySortingdir);
%         if ~isempty(dir(fullfile(mySortingdir, sprintf('Record Node*'))))
%             tempfolder = (dir(fullfile(mySortingdir, sprintf('Record Node*'))));
%             mySortingdir = fullfile(mySortingdir,tempfolder.name);
%         end
%     else
%     % more recording files than behavior files - send both
%         for i = 1:size(myfolders,1)
%             ephysdir = myfolders(i).name;
%             ephysdir = fullfile(root,ephysdir);
%             if ~isempty(dir(fullfile(ephysdir, sprintf('Record Node*'))))
%                 tempfolder = (dir(fullfile(ephysdir, sprintf('Record Node*'))));
%                 ephysdir = fullfile(ephysdir,tempfolder.name);
%             end
%             mySortingdir(i,:) = ephysdir;
%         end
%     end
end

end