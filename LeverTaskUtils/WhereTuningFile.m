function [TuningFile] = WhereTuningFile(FilePaths,BehaviorFile)
TuningFile = [regexprep(BehaviorFile,'_r','_o'),'.mat'];
if exist(fullfile(FilePaths,TuningFile),'file')
    TuningFile = fullfile(FilePaths,TuningFile);
else
    TuningFile = [];
end
end