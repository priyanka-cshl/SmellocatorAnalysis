function [TuningFile] = WhereTuningFile(FilePaths,BehaviorFile)
TuningFile = [regexprep(BehaviorFile,'_r','_o'),'.mat'];
if exist(fullfile(FilePaths,TuningFile),'file')
    
    % check if there's more than one tuning file - mid session Passive replays
    FileList = dir([fullfile(FilePaths,TuningFile(1:end-5)),'*']);
    if numel(FileList)>1
        % is there also more than one behavior session?
        BehaviorList = dir([fullfile(FilePaths,BehaviorFile(1:end-1)),'*']);
        if numel(BehaviorList)>1
            % likely a mid-session passive replay session
            clear TuningFile
            for n = 1:numel(FileList)
                TuningFile(n,:) = fullfile(FilePaths,FileList(n).name);
            end
        end
    else
        TuningFile = fullfile(FilePaths,TuningFile);
    end
    
else
    TuningFile = [];
end
end