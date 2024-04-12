MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S1/S1_20230327_r0_processed.mat';
MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q4/Q4_20221109_r0_processed.mat';
MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q8/Q8_20221204_r0_processed.mat';

[TrialAligned, TrialInfo, ReplayAligned, ReplayInfo, TuningAligned, TuningInfo, AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);


% get all sniffs of a particular stimulus state
% State = [odorstate location CL/OL/Tuning perturbed]
% odorstate = -1 (no air), 0 (Air), 1,2,3 (odors)
% location  = -100 to 100
% arrayfun(@(x) find(x{1}(:,2)==-1), TrialAligned.Sniffs, 'UniformOutput', false)

%%
odorstate   = 1; % IAA
whichUnit = 72;
thisUnitSpikes = AllUnits.Spikes{whichUnit};

%% for closed loop
% get all unpertubed trials
whichtrials = find(cellfun(@isempty, TrialInfo.Perturbation(:,1))); % vanilla closed-loop
% also include trials marked as OL-Template - they are also just control
% trials
whichtrials = sort([whichtrials; find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'))]); % also used the template trials (no-perturbations)

% extract all sniffs that happen during trials of the selected odor
AllSniffs = [];
for i = 1:numel(whichtrials)
    thisTrial       = whichtrials(i); % Trial ID
    thisTrialStart  = TrialAligned.RefTS(thisTrial,1);
    
    whichsniffs     = find(TrialAligned.Sniffs{thisTrial}(:,2)==odorstate);
    
    SniffTimestamps = TrialAligned.Sniffs{thisTrial}(whichsniffs,:);
    SniffTimestamps(:,5:11) = SniffTimestamps(:,5:11) + thisTrialStart;
    AllSniffs       = [AllSniffs; SniffTimestamps];
end

% add a column for current sniff/inhalation duration
AllSniffs(:,15) = AllSniffs(:,8) - AllSniffs(:,7);
AllSniffs = sortrows(AllSniffs, [4 15]); % snifftype and duration

% % sort sniffs by odor location (col 13) and/or duration
% AllSniffs = sortrows(AllSniffs, [13 15]);

% get sniff aligned Spikes for the user selected unit

figure;

subplot(1,3,1);
hold on

% for sniff periods
PlotExhalationPeriods(AllSniffs);

% for spikes
[SpikesPSTH, YMax] = SniffAlignedSpikesPlotter(AllSniffs,thisUnitSpikes);

set(gca,'XLim', [-0.1 1]);

%%
% for replays
% get replay trials
whichtrials = 1:size(ReplayInfo.TrialID,1);

% extract all sniffs that happen during trials of the selected odor
AllSniffs = [];
for i = 1:numel(whichtrials)
    thisTrial       = whichtrials(i); % Trial ID
    thisTrialStart  = ReplayAligned.RefTS(thisTrial,1);
    
    whichsniffs     = find(ReplayAligned.Sniffs{thisTrial}(:,2)==odorstate);
    
    SniffTimestamps = ReplayAligned.Sniffs{thisTrial}(whichsniffs,:);
    SniffTimestamps(:,5:11) = SniffTimestamps(:,5:11) + thisTrialStart;
    AllSniffs       = [AllSniffs; SniffTimestamps];
end

% add a column for current sniff/inhalation duration
AllSniffs(:,15) = AllSniffs(:,8) - AllSniffs(:,7);
AllSniffs = sortrows(AllSniffs, [15]);

% % sort sniffs by odor location (col 13) and/or duration
% AllSniffs = sortrows(AllSniffs, [13 15]);

subplot(1,3,2);
hold on
% plot the sniff markers
SniffIndexes = (1:size(AllSniffs,1))';
SniffIndexes = reshape(repmat(SniffIndexes',2,1),[],1);
SniffIndexes = circshift(SniffIndexes,-1);
SniffIndexes(end) = SniffIndexes(end-1) + 1;

for col = 5:10
    SniffBounds = AllSniffs(:,col)-AllSniffs(:,7);
    plot(reshape(repmat(SniffBounds',2,1),[],1),SniffIndexes,'r');
end

for s = 1:size(AllSniffs,1)
    t = AllSniffs(s,[5 7 11]); % [prevsniffstart thisniffstart nextsniffend]
    whichSpikes = intersect(find(AllUnits.Spikes{whichUnit}>=t(1)), find(AllUnits.Spikes{whichUnit}<=t(3)));
    thisSniffSpikes = AllUnits.Spikes{whichUnit}(whichSpikes) - t(2);
    if ~isempty(thisSniffSpikes)
        plot(thisSniffSpikes,s*(ones(numel(thisSniffSpikes))),'*k');
    end
end
set(gca,'YLim', [-1 YMax+1]);

% for tuning
% get tuning trials
whichtrials = 1:size(TuningInfo.TrialID,1);

% extract all sniffs that happen during trials of the selected odor
AllSniffs = [];
for i = 1:numel(whichtrials)
    thisTrial       = whichtrials(i); % Trial ID
    thisTrialStart  = TuningAligned.RefTS(thisTrial,1);
    
    whichsniffs     = find(TuningAligned.Sniffs{thisTrial}(:,2)==odorstate);
    
    SniffTimestamps = TuningAligned.Sniffs{thisTrial}(whichsniffs,:);
    SniffTimestamps(:,5:11) = SniffTimestamps(:,5:11) + thisTrialStart;
    AllSniffs       = [AllSniffs; SniffTimestamps];
end

% add a column for current sniff/inhalation duration
AllSniffs(:,15) = AllSniffs(:,8) - AllSniffs(:,7);
AllSniffs = sortrows(AllSniffs, [15]);

% % sort sniffs by odor location (col 13) and/or duration
% AllSniffs = sortrows(AllSniffs, [13 15]);

if ~isempty(AllSniffs)
    subplot(1,3,3);
    hold on
    % plot the sniff markers
    SniffIndexes = (1:size(AllSniffs,1))';
    SniffIndexes = reshape(repmat(SniffIndexes',2,1),[],1);
    SniffIndexes = circshift(SniffIndexes,-1);
    SniffIndexes(end) = SniffIndexes(end-1) + 1;
    
    for col = 5:10
        SniffBounds = AllSniffs(:,col)-AllSniffs(:,7);
        plot(reshape(repmat(SniffBounds',2,1),[],1),SniffIndexes,'r');
    end
    
    for s = 1:size(AllSniffs,1)
        t = AllSniffs(s,[5 7 11]); % [prevsniffstart thisniffstart nextsniffend]
        whichSpikes = intersect(find(AllUnits.Spikes{whichUnit}>=t(1)), find(AllUnits.Spikes{whichUnit}<=t(3)));
        thisSniffSpikes = AllUnits.Spikes{whichUnit}(whichSpikes) - t(2);
        if ~isempty(thisSniffSpikes)
            plot(thisSniffSpikes,s*(ones(numel(thisSniffSpikes))),'*k');
        end
    end
    set(gca,'YLim', [-1 YMax+1]);
end

function [] = PlotExhalationPeriods(AllSniffs,whichcols)

% plot the sniff markers
SniffIndexes = (1:size(AllSniffs,1))';
SniffIndexes = reshape(repmat(SniffIndexes',2,1),[],1);
SniffIndexes = circshift(SniffIndexes,-1);
SniffIndexes(end) = SniffIndexes(end-1) + 1;

if nargin<2
    whichcols = 0;
end

% whichcols = [-1 0 1]
for ncols = 1:numel(whichcols)
    mycol = 8 + 2*whichcols(ncols); % [6 8 10]
    exh_markers = [];
    SniffBounds = AllSniffs(:,mycol)-AllSniffs(:,7);
    exh_markers = [reshape(repmat(SniffBounds',2,1),[],1) SniffIndexes];
    SniffBounds = AllSniffs(:,mycol+1)-AllSniffs(:,7);
    exh_markers = [exh_markers;
                flipud([reshape(repmat(SniffBounds',2,1),[],1) SniffIndexes]) ];
    patch(exh_markers(:,1),exh_markers(:,2),Plot_Colors('TZ'),...
        'FaceAlpha',0.3, ...
        'EdgeAlpha',0.3, ...
        'edgecolor',Plot_Colors('TZ'));
end

end

function [SpikesPSTH, YMax] = SniffAlignedSpikesPlotter(AllSniffs,thisUnitSpikes,YMax)
SpikesPlot = []; SpikesPSTH = [];
for s = 1:size(AllSniffs,1)
    t = AllSniffs(s,[5 7 11 9]); % [prevsniffstart thisniffstart nextsniffend thissniffend]
    whichSpikes = intersect(find(thisUnitSpikes>=t(1)), find(thisUnitSpikes<=t(3)));
    thisSniffSpikes = thisUnitSpikes(whichSpikes) - t(2);
    SpikesPlot = vertcat(SpikesPSTH, [thisSniffSpikes s*ones(numel(thisSniffSpikes),1)]);
    % for psth - might want to ignore some periods
    thisSniffSpikes(:,thisSniffSpikes>t(4)) = [];
    SpikesPSTH = vertcat(SpikesPSTH, [thisSniffSpikes s*ones(numel(thisSniffSpikes),1)]);
end

if ~isempty(SpikesPlot)
    plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize',0.5);
end

if nargin<3
    YMax = s;
end
set(gca,'YLim', [-1 YMax+1]);

end