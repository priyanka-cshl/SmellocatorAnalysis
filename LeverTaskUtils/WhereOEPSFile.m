function [myephysdir] = WhereOEPSFile(BehaviorFile,BehaviorFolder)
foo = strsplit(BehaviorFile,'_');
mousename = char(foo(1));
date = char(foo(2));
datetoken = [date(1:4),'-',date(5:6),'-',date(7:8)];
    
[Paths] = WhichComputer();
if exist(fullfile(Paths.Grid.Ephys{1},mousename),'dir')
    root = fullfile(Paths.Grid.Ephys{1},mousename);
elseif exist(fullfile([Paths.Grid.Ephys{2},mousename(1:end-1)],mousename))
    root = fullfile([Paths.Grid.Ephys{2},mousename(1:end-1)],mousename);
elseif exist(fullfile(Paths.Grid.Ephys{3},mousename))
    root = fullfile(Paths.Grid.Ephys{3},mousename);
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
    % same # of behavior and recording files
    if numel(dir(fullfile(BehaviorFolder,['*',date,'*r*.mat']))) == size(myfolders,1)
        sessionname = char(foo(3));
        session_num = str2num(sessionname(2));
        myephysdir = myfolders(session_num+1).name;
        myephysdir = fullfile(root,myephysdir);
        if ~isempty(dir(fullfile(myephysdir, sprintf('Record Node*'))))
            tempfolder = (dir(fullfile(myephysdir, sprintf('Record Node*'))));
            myephysdir = fullfile(myephysdir,tempfolder.name);
        end
    else
    % more recording files than behavior files - send both
        for i = 1:size(myfolders,1)
            ephysdir = myfolders(i).name;
            ephysdir = fullfile(root,ephysdir);
            if ~isempty(dir(fullfile(ephysdir, sprintf('Record Node*'))))
                tempfolder = (dir(fullfile(ephysdir, sprintf('Record Node*'))));
                ephysdir = fullfile(ephysdir,tempfolder.name);
            end
            myephysdir(i,:) = ephysdir;
        end
    end
end

end