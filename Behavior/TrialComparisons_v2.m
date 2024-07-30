%% script to compare movement patterns across trials and across target zones

%% Paths
%Mouse = 'Q5'; Date = '20221103';
DataRoot    = '/home/priyanka/Dropbox/Smellocator_Behavior_2/';
Session = dir(fullfile(DataRoot,Mouse,['*',Date,'*']));
SessionPath = fullfile(DataRoot,Mouse,Session.name);

% Load the relevant variables
load(SessionPath, 'Traces', 'TrialInfo', 'TargetZones', 'SniffTS',...
               'startoffset', 'SampleRate');

%% 
AllTargets  = unique(TrialInfo.TargetZoneType);
AllOdors    = unique(TrialInfo.Odor);
TargetLims  = 1:0.25:3.75;
Pretrial    = 0.2; % in seconds

figure;
[ha, pos] = tight_subplot(numel(AllTargets),numel(AllOdors),[.01 .01],[.01 .01],[.01 .01]);
xlength = [];
for i = 1:numel(AllTargets)
    for j = 1:numel(AllOdors)
        % select subplot
        axes(ha((3*i)-3+j));
        hold on
        plot_ts     = 0;
        
        % all sniffs of this target zone and this odor
        f = intersect( find(TrialInfo.TargetZoneType == AllTargets(i)) , ...
                    find(TrialInfo.Odor == AllOdors(j)) );
                
        for k = 1:numel(f) % every valid trial
            trialID = f(k);
            idx(1) = find(diff((Traces.TrialState{trialID}))>0,1,'first'); % trial start
            idx(2) = find(diff((Traces.TrialState{trialID}))<0,1,'last'); % trial stop
            
            idx(1) = idx(1) - (Pretrial*SampleRate); % take some samples before trial start
            
            TraceTimeStamps = Traces.Timestamps{trialID}(idx(1):idx(2));
            if ~isempty(SniffTS)
                thisTrialSniffs = SniffTS(find(SniffTS(:,5)==trialID),:);
                % remove sniffs that happened outside the trace window
                thisTrialSniffs(thisTrialSniffs(:,1)<TraceTimeStamps(1),:) = [];
                thisTrialSniffs(thisTrialSniffs(:,1)>TraceTimeStamps(end),:) = [];
            else
                % detect sniff timestamps per trial
                Thermistor      = Traces.Sniffs{trialID}(idx(1):idx(2));
                thisTrialSniffs = ProcessThermistorData([TraceTimeStamps Thermistor]);
            end
            LeverTrace = [];
            LeverTrace(1,:) = Traces.Lever{trialID}(idx(1):idx(2));
            LeverTrace(2,:) = LeverTrace(1,:);
            
            for n = 1:size(thisTrialSniffs,1)
                inhalation_indices = intersect( ...
                    find(TraceTimeStamps>=thisTrialSniffs(n,1)), ...
                    find(TraceTimeStamps<thisTrialSniffs(n,2)) );
                LeverTrace(2,inhalation_indices) = NaN;
            end
            
            PlotTimeStamps = TraceTimeStamps - TraceTimeStamps(1) + plot_ts;
            plot_ts = PlotTimeStamps(end) + Pretrial;
            
            PlotTrialPeriod([(PlotTimeStamps(1) + Pretrial) PlotTimeStamps(end) TargetLims(TrialInfo.TargetZoneType(trialID))]);

            plot(PlotTimeStamps,LeverTrace(1,:),'r','LineWidth',1.5 - ~TrialInfo.Success(trialID));
            plot(PlotTimeStamps,LeverTrace(2,:),'k','LineWidth',1.5 - ~TrialInfo.Success(trialID));
            
            xlength((3*i)-3+j) = PlotTimeStamps(end) + Pretrial;
        end
        set(gca, 'YLim', [0 5.5]);
    end
end

% rescale time axis of all plots
for i = 1:numel(xlength)
    set(ha(i),'XLim',[0 max(xlength)]);
end