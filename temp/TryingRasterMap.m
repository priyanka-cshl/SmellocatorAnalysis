[WhichSession, SessionPath] = uigetfile(...
                                fullfile('O3/O3_20210922_r0_processed.mat'),...
                                'Select Behavior or Recording Session');
WhereSession = fullfile(SessionPath,WhichSession);

% Load the relevant variables
load(WhereSession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
               'startoffset', 'errorflags', 'SampleRate', ...
               'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');
           
if ~isempty(TTLs)
    SessionLength = num2str(10*ceil(TTLs.Trial(end,2)/10));
    NumUnits = size(SingleUnits,2);
else
    SessionLength = TrialInfo.SessionTimestamps(end,2);
    NumUnits = NaN;
end

if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    u = unique(TrialInfo.Perturbation(x));
    PerturbationList = u{1};
    for y = 2:size(u,1)
        PerturbationList = [PerturbationList,'; ',u{y}];
    end
else
    x = [];
    PerturbationList = '';
end

%% Concatenate traces
whichTrials = 1:length(TrialInfo.TrialID);
traceOverlap = SampleRate*startoffset;
%whichTraces = fieldnames(Traces);
whichTraces{1} = 'Lever'; whichTraces{2} = 'Motor'; whichTraces{3} = 'Sniffs';
whichTraces{4} = 'Trial'; whichTraces{5} = 'Rewards'; whichTraces{6} = 'Timestamps';

for j = 1:size(whichTraces,2)
    temp = cellfun(@(x) ...
        x(1:end-traceOverlap), Traces.(whichTraces{j})(whichTrials), ...
        'UniformOutput', false);
    TracesOut.(whichTraces{j}) = {[cell2mat(temp(:)); ...
        Traces.(whichTraces{j}){whichTrials(end)}(end-traceOverlap+1:end,1)]};
end

% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);

if ~isempty(TTLs)
    TrialStart_Ephys = TTLs.Trial(1,2);
    % factor to convert all behavior timestamps to match Ephys
    TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;
else
    TrialStart_Ephys = 0;
    TimestampAdjuster = 0;
end

TracesOut.Timestamps{1} = Timestamps + TimestampAdjuster;

%% Spike data
TS = 1000*[447 497]; % 50 seconds = 50,000 bins
timeBins = TS(1):TS(2);
myRaster = zeros(NumUnits,numel(timeBins));
TS = TS/1000;

figure; hold on;

% plot Odor boxes
BoxTag = 'odorodor';
for i = 1:3
    handles.([BoxTag,num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.([BoxTag,num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TTLs.(['Odor',num2str(i)])(:,1:2)';
    ValveTS(:,1:(find(ValveTS(1,:)<=TS(1),1,'last'))) = [];
    ValveTS(:,(find(ValveTS(2,:)>=TS(2),1,'first')):end) = [];
    ValveTS = ValveTS - TS(1);
    handles.([BoxTag,num2str(i),'Plot']).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat((10+NumUnits)*[0 1 1 0]',size(ValveTS,2),1)];
    handles.([BoxTag,num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end

for n = 1:NumUnits
    % align to the specified event
    f = find(SingleUnits(n).spikes>TS(2),1,'first');
    if ~isempty(f)
        f = f-1;
        thisUnitSpikes = SingleUnits(n).spikes(1:f)' - TS(1);
    else
        thisUnitSpikes = SingleUnits(n).spikes' - TS(1);
    end
    PlotRaster(thisUnitSpikes(thisUnitSpikes>0),n,[0 0 0],0.8);
    % convert spike times to milliseconds and floor values
	thisUnitSpikes = floor(1000*thisUnitSpikes);
    % remove NaNs
    thisUnitSpikes(isnan(thisUnitSpikes)) = [];
	% Make raster
    [C,~,ic] = unique(thisUnitSpikes);
    bin_counts = accumarray(ic,1);
    if ~isempty(C)
        % ignore any -ve time bins
        bin_counts((C<=0),:) = [];
        C(C<=0) = [];
        myRaster(n,C) = bin_counts;
    end
end

%% run rastermap
ops.nCall = [30 NumUnits];
ops.iPC = 1:NumUnits;
[isort1, isort2, Sm] = mapTmap(myRaster, ops);

%% replot
figure; 
hold on;

% plot Odor boxes
BoxTag = 'odor';
for i = 1:3
    handles.([BoxTag,num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.([BoxTag,num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TTLs.(['Odor',num2str(i)])(:,1:2)';
    ValveTS(:,1:(find(ValveTS(1,:)<=TS(1),1,'last'))) = [];
    ValveTS(:,(find(ValveTS(2,:)>=TS(2),1,'first')):end) = [];
    ValveTS = ValveTS - TS(1);
    handles.([BoxTag,num2str(i),'Plot']).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat((10+NumUnits)*[0 1 1 0]',size(ValveTS,2),1)];
    handles.([BoxTag,num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end

for m = 1:NumUnits
    n = isort1(m);
    % align to the specified event
    f = find(SingleUnits(n).spikes>TS(2),1,'first');
    if ~isempty(f)
        f = f-1;
        thisUnitSpikes = SingleUnits(n).spikes(1:f)' - TS(1);
    else
        thisUnitSpikes = SingleUnits(n).spikes' - TS(1);
    end
    PlotRaster(thisUnitSpikes(thisUnitSpikes>0),m,[0 0 0],0.8);
end


%% replot
figure; 
hold on;
TS = [1027 1077];

% plot Odor boxes
BoxTag = 'Rodor';
for i = 1:3
    handles.([BoxTag,num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.([BoxTag,num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TTLs.(['Odor',num2str(i)])(:,1:2)';
    ValveTS(:,1:(find(ValveTS(1,:)<=TS(1),1,'last'))) = [];
    ValveTS(:,(find(ValveTS(2,:)>=TS(2),1,'first')):end) = [];
    ValveTS = ValveTS - TS(1);
    handles.([BoxTag,num2str(i),'Plot']).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat((10+NumUnits)*[0 1 1 0]',size(ValveTS,2),1)];
    handles.([BoxTag,num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end

for m = 1:NumUnits
    n = isort1(m);
    % align to the specified event
    f = find(SingleUnits(n).spikes>TS(2),1,'first');
    if ~isempty(f)
        f = f-1;
        thisUnitSpikes = SingleUnits(n).spikes(1:f)' - TS(1);
    else
        thisUnitSpikes = SingleUnits(n).spikes' - TS(1);
    end
    PlotRaster(thisUnitSpikes(thisUnitSpikes>0),m,[0 0 0],0.8);
end
