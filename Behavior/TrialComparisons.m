%% script to compare movement patterns across trials and across target zones

%% Paths
DataRoot    = '/home/priyanka/Dropbox/Smellocator_Behavior/Q3';
SessionName = 'Q3_20221005_r0_processed.mat';
SessionPath = fullfile(DataRoot,SessionName);

% Load the relevant variables
load(SessionPath, 'Traces', 'TrialInfo', 'TargetZones', 'SniffTS',...
               'startoffset', 'SampleRate');

%% check if sniff time stamps have been processed
if ~exist('SniffTS','var')
    [Paths] = WhichComputer();
    [~,Mouse] = fileparts(DataRoot);
    RawData = fullfile(Paths.Grid.Behavior,Mouse,strrep(SessionName,'_processed',''));
    [SniffTS] = ReadThermistorData(RawData); % in behavior timestamps    
end

%% 
AllTargets  = unique(TrialInfo.TargetZoneType);
AllOdors    = unique(TrialInfo.Odor);
for i = 1:numel(AllTargets)
    for j = 1:numel(AllOdors)
        f = intersect( find(TrialInfo.TargetZoneType == AllTargets(i)) , ...
                    find(TrialInfo.Odor == AllOdors(j)) );
                
        for k = 1:numel(f) % every valid trial
            trialID = f(k);
            idx(1) = find(diff((Traces.TrialState{trialID}))>0,1,'first');
            idx(2) = find(diff((Traces.TrialState{trialID}))<0,1,'first');
            
            idx(1) = idx(1) - (1000/SampleRate)*100;
            
            thisTrialSniffs = intersect( ...
                find(SniffTS(:,1)>=Traces.Timestamps{trialID}(idx(1))) , ...
                find(SniffTS(:,1)< Traces.Timestamps{trialID}(idx(2))) );
            thisTrialSniffs = SniffTS(thisTrialSniffs,:);
            
            LeverTrace(1,:) = Traces.Lever{trialID}(idx(1):idx(2));
            LeverTrace(2,:) = LeverTrace(1,:);
            TraceTimeStamps = Traces.Timestamps{trialID}(idx(1):idx(2));
            for n = 1:2:size(thisTrialSniffs,1)
                sniff_indices = intersect( ...
                    find(TraceTimeStamps>=thisTrialSniffs(n,1)), ...
                    find(TraceTimeStamps<thisTrialSniffs(n,3)) );
                LeverTrace(2,sniff_indices) = NaN;
            end
        end
    end
end