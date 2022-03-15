%% for O3
SessionPath = 'O3/O3_20210927_r0_processed.mat';
% Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
AlignTo = [2 2 2]; 
ChosenUnits = [18 25 53];
HaltUnits = [63 15 10]; % 31 17]
HaltUnits = [31 17 18];
FRmax = [0 50];

%% for O3
SessionPath = 'O3/O3_20210929_r0_processed.mat';
% Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
AlignTo = [2 2 2]; 
ChosenUnits = [18 25 53];
HaltUnits = [16 74 15];
FRmax = [0 50];

% %% for PCX4
% SessionPath = 'PCX4/PCX4_20210713_r0_processed.mat';
% % Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
% AlignTo = [1 1 1]; 
% ChosenUnits = [7 15 43];
% HaltUnits = [7 17 18];
% FRmax = [0 20];

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

%% Plotting example responses
nCols = numel(HaltUnits);
figure('Position',[0 0 1000 600]);

for x = 1:nCols
    whichUnit = HaltUnits(x);
    
    trialsdone = 0;
    % control trials
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
    end
    
    % halt trials
    for whichTZ = 1:4:12
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

