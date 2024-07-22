function [myOEPSdir] = WhereRawOEPS(BehaviorFile,BehaviorFolder)
foo = strsplit(BehaviorFile,'_');
mousename = char(foo(1));
date = char(foo(2));
datetoken = [date(1:4),'-',date(5:6),'-',date(7:8)];
    
[Paths] = WhichComputer();

if exist(fullfile(Paths.Grid.Ephys{3},mousename),'dir')
    root = fullfile(Paths.Grid.Ephys{3},mousename);
end

myfolders = dir ([root,'/',datetoken,'*']);

if size(myfolders,1) == 0
    myOEPSdir = [];
    disp('no oeps file found');
elseif size(myfolders,1) == 1 % only one recording session found
    myOEPSdir = myfolders.name;
    myOEPSdir = fullfile(root,myOEPSdir);
else

end