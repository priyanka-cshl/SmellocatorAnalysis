%% for O3
SessionPath = 'O3/O3_20210927_r0_processed.mat';
% Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
AlignTo = [2 2 2]; 
ChosenUnits = [18 25 53];
HaltUnits = [63 15 10]; % 31 17]
FRmax = [0 30];
 
% %% for PCX5
% SessionPath = 'PCX5/PCX5_20210713_r0_processed.mat';
% % Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
% AlignTo = [1 5 4]; 
% ChosenUnits = [30 34 43];
% FRmax = [0 50];

%% for PCX4
SessionPath = 'PCX4/PCX4_20210713_r0_processed.mat';
% Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
AlignTo = [1 1 1]; 
ChosenUnits = [7 15 43];
FRmax = [0 50];

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');
            
[~, ~, TrialInfo] = LoadProcessedDataSession(MySession); % to get target zone entry time

%% Trial Aligned spikeTimes
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

%%
whichOdor = 1; 
for tz = 1:12
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
                    find(TrialInfo.Odor==whichOdor));
    whichTrials = intersect(whichTrials,...
                    find(TrialInfo.TargetZoneType==tz));
    OdorONPeriod(tz,1) = median(TrialInfo.OdorStart(whichTrials,1));
    foo = cellfun(@(x) median(x(491:500)), Traces.Motor(whichTrials),'UniformOutput', false);
    LocationTrialStart(tz,1) = round(median(cell2mat(foo)));
end

%% Plot
nCols = numel(ChosenUnits);
figure('Position',[0 0 1000 600]);
for x = 1:nCols
    whichUnit = ChosenUnits(x);
    
    switch AlignTo(x)
        case {1,2,6}
            myXlim = [-1.2 5];
        case {3,4}
            myXlim = [-5.2 1];
        case 5
            myXlim = [-2.2 5];
    end

    trialsdone = 0;
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        hold on
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter(whichUnit, whichOdor, whichTZ, AlignedSpikes, Events, TrialInfo, AlignTo(x), trialsdone, trialboxcolor);
        
        set(gca, 'XLim', myXlim, 'TickDir', 'out');
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        hold on
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);

        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        
        trialsdone = trialsdone + nTrials;
    end    
end

%% Plotting example responses
nCols = numel(HaltUnits);
figure('Position',[0 0 1000 600]);
for x = 1:nCols
    whichUnit = HaltUnits(x);
    
    trialsdone = 0;
    for whichTZ = 1:4:12
        
        trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        hold on
        
        % first plot control trials
        % Spikes
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter(whichUnit, whichOdor, whichTZ, AlignedSpikes, Events, TrialInfo, 2, trialsdone, trialboxcolor);
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        hold on
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        % halt trials
        trialboxcolor = GetLocationColor(-LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter(whichUnit, whichOdor, -whichTZ, AlignedSpikes, Events, TrialInfo, 6, trialsdone, trialboxcolor);
        
        set(gca, 'XLim', myXlim, 'TickDir', 'out');
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        
        
    end
end

%% Behavior plotting
for i = 1:12
    whichTZ = i;
    
    % control trials
    controlTrials = intersect(find(TrialInfo.TargetZoneType==whichTZ),...
                              find(~strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')));
    controlTrials(:,3) = Events(controlTrials,5);

    for j = 1:numel(controlTrials(:,1))
        Traces.HaltFlip(controlTrials(j,1)) = {0*Traces.Trial{controlTrials(j,1)}};
        Traces.TargetZone(controlTrials(j,1)) = {TrialInfo.TargetZoneType(controlTrials(j,1)) + ...
            0*Traces.Trial{controlTrials(j,1)}};
    end
        
    % perturbation trials
    perturbationTrials = intersect(find(TrialInfo.TargetZoneType==whichTZ),...
                              find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')));
    perturbationTrials(:,3) = Events(perturbationTrials,5);
    
    for j = 1:numel(perturbationTrials(:,1))
        foo = 0*Traces.Trial{perturbationTrials(j,1)};
        idx = TrialInfo.Perturbation{perturbationTrials(j,1),2};
        idx(:,1:2) = idx(:,1:2) + SampleRate*startoffset;
        foo(idx(1):idx(2)) = -idx(3);
        Traces.HaltFlip(perturbationTrials(j,1)) = {foo};
        Traces.TargetZone(perturbationTrials(j,1)) = ...
                {TrialInfo.TargetZoneType(perturbationTrials(j,1)) + ...
                0*foo};
    end
end

B1 = figure;
figure(B1);
[TracesOut] = ConcatenateTraces(Traces, 1:size(Traces.Lever,2), SampleRate*startoffset);

timestamps = (1:size(TracesOut.Lever{1},1))'/SampleRate;
TZ = [TracesOut.TargetZone{1} TracesOut.TargetZone{1}];
TZ = [TargetZones(TZ(:,1),1) TargetZones(TZ(:,1),3)];
Trial = TracesOut.Trial{1};
%Trial(Trial~=whichOdor) = 0;
PlotBehaviorPerturbations(timestamps,TracesOut.Lever{1},...
    TracesOut.Sniffs{1},TracesOut.Licks{1},TracesOut.Rewards{1},...
    Trial,...
    TZ, ...
    TracesOut.HaltFlip{1},5);

% make odor ON time available on the plot
traceoffset = ((TrialInfo.TraceIndices(1,1))-1)/SampleRate;
plot(-traceoffset+(TrialInfo.OdorStart(:,1)+TrialInfo.SessionTimestamps(:,1)),-0.5,'.k');      
set(gcf,'Position',[0 100 1000 300]);

%% for O3
% set(gca,'XLim',[725 736],'YLim',[-1 8],'TickDir','out')
% set(gca,'XLim',[961.9 972.9],'YLim',[-1 8],'TickDir','out')