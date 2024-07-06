function [SingleUnits] = ReSaveSingleUnits(SessionPath)

load(SessionPath, 'FileLocations', 'TTLs', 'TuningTTLs', 'TrialInfo');
if isfield(FileLocations,'Spikes')
    SingleUnits = GetSingleUnits(FileLocations.Spikes);
    [SingleUnits] = Spikes2Trials(SingleUnits, TTLs.Trial(1:size(TrialInfo.TrialID,2),:), TuningTTLs);
    save(SessionPath, 'SingleUnits', '-append');
else
    disp('warning: this session has not been processed before for SingleUnits');
    keyboard;
end

end
