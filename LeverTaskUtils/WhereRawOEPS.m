function [myOEPSdir] = WhereRawOEPS(BehaviorFile,BehaviorFolder)
foo = strsplit(BehaviorFile,'_');
mousename = char(foo(1));
date = char(foo(2));
datetoken = [date(1:4),'-',date(5:6),'-',date(7:8)];
    
[Paths] = WhichComputer();

rootdir = [];
if exist(fullfile(Paths.Grid.Ephys{3},mousename),'dir')
    rootdir = fullfile(Paths.Grid.Ephys{3},mousename);
elseif exist(fullfile([Paths.Grid.Ephys{2},mousename(1)],mousename),'dir')
    rootdir = fullfile([Paths.Grid.Ephys{2},mousename(1)],mousename);
end

myfolders = dir ([rootdir,'/',datetoken,'*']);

if size(myfolders,1) == 0
    myOEPSdir = [];
    disp('no oeps file found');
elseif size(myfolders,1) == 1 % only one recording session found
    myOEPSdir = myfolders.name;
    myOEPSdir = fullfile(rootdir,myOEPSdir);
else

end