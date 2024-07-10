%% script to compare movement patterns across trials and across target zones

%% Paths
%DataRoot    = '/home/priyanka/Dropbox/Smellocator_Behavior/Q3';
DataRoot    = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Ryan_Behavior/Q4';
SessionName = 'Q4_20221103_r0_processed.mat';
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
TargetLims  = 1:0.25:3.75;
for i = 1:numel(AllTargets)
    for j = 1:numel(AllOdors)
        subplot(numel(AllTargets),numel(AllOdors),(3*i)-3+j);
        hold on
        plot_ts     = 0;

        f = intersect( find(TrialInfo.TargetZoneType == AllTargets(i)) , ...
                    find(TrialInfo.Odor == AllOdors(j)) );
                
        for k = 1:numel(f) % every valid trial
            trialID = f(k);
            idx(1) = find(diff((Traces.TrialState{trialID}))>0,1,'first');
            idx(2) = find(diff((Traces.TrialState{trialID}))<0,1,'first');
            
            idx(1) = idx(1) - (0.1*SampleRate);
            
            thisTrialSniffs = intersect( ...
                find(SniffTS(:,1)>=Traces.Timestamps{trialID}(idx(1))) , ...
                find(SniffTS(:,1)< Traces.Timestamps{trialID}(idx(2))) );
            thisTrialSniffs = SniffTS(thisTrialSniffs,:);
            
            LeverTrace = [];
            LeverTrace(1,:) = Traces.Lever{trialID}(idx(1):idx(2));
            LeverTrace(2,:) = LeverTrace(1,:);
            TraceTimeStamps = Traces.Timestamps{trialID}(idx(1):idx(2));
            for n = 1:2:size(thisTrialSniffs,1)
                sniff_indices = intersect( ...
                    find(TraceTimeStamps>=thisTrialSniffs(n,1)), ...
                    find(TraceTimeStamps<thisTrialSniffs(n,3)) );
                LeverTrace(2,sniff_indices) = NaN;
            end
            
            PlotTimeStamps = TraceTimeStamps - TraceTimeStamps(1) + plot_ts;
            plot_ts = PlotTimeStamps(end) + 0.1;
            
            PlotTrialPeriod([(PlotTimeStamps(1) + 0.1) PlotTimeStamps(end) TargetLims(TrialInfo.TargetZoneType(trialID))]);

            plot(PlotTimeStamps,LeverTrace(1,:),'k','LineWidth',1+TrialInfo.Success(trialID));
            plot(PlotTimeStamps,LeverTrace(2,:),'r','LineWidth',1+TrialInfo.Success(trialID));

        end
        set(gca, 'YLim', [0 5]);
    end
end