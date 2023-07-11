function [myimagingdir] = WhereWidefieldData(BehaviorFile,AnimalName)
[Paths] = WhichComputer();
foo = regexprep(BehaviorFile,[AnimalName,'_'],'');
myimagingdir = fullfile(Paths.Widefield.Processed,AnimalName,foo);

if ~exist(myimagingdir)
    myimagingdir = [];
    disp('no imging directory found');
end

end