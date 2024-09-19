
function [TracesOut,SingleUnits] = LoadLeverTaskSession_wdw(WhereSession,savemode)

if nargin < 2
    savemode = 0;
end

% load the processed behavior/ephys files (from PreprocessSmellocatorData.m)
% [WhichSession, SessionPath] = uigetfile(...
%                                 fullfile('Q4/Q4_20221109_r0_processed.mat'),...
%                                 'Select Behavior or Recording Session');
% WhereSession = fullfile(SessionPath,WhichSession);

% WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q4/Q4_20221109_r0_processed.mat';
% WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/O3/O3_20210929_r0_processed.mat';
% WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q3/Q3_20221019_r0_processed.mat';
% WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S12/S12_20230731_r0_processed.mat';
% WhereSession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q9/Q9_20221116_r0_processed.mat';

%% Load the relevant variables
load(WhereSession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits', ...
    'TimestampAdjust', 'FileLocations');

%% flag perturbation trials
if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    perturbedTrials = x(~strcmp(TrialInfo.Perturbation(x,1),'OL-Template'));
else
    perturbedTrials = [];
end

%% Change the trial vector to include the OdorON duration before trial start
for i = 1:size(TrialInfo.TrialID,2)
    OdorTrace = Traces.Trial{i};
    idx2 = find(OdorTrace((SampleRate*startoffset):end)>0,1,'first') + (SampleRate*startoffset) - 1;
    idx1 = idx2 + round(TrialInfo.OdorStart(i,1)*SampleRate);

    if idx1<1
        %keyboard;

        % for the current trial
        OdorTrace(OdorTrace>0) = 1;
        OdorTrace(1:idx2,1) = 1;
        Traces.Odor{i} = OdorTrace*TrialInfo.Odor(i,1);

        % previous trial
        prevOdorTrace = Traces.Odor{i-1};
        prevOdorTrace(end+idx1:end) = TrialInfo.Odor(i,1);
        Traces.Odor{i-1} = prevOdorTrace;
    else
        OdorTrace(OdorTrace>0) = 1;
        OdorTrace(idx1:idx2,1) = 1;
        Traces.Odor{i} = OdorTrace*TrialInfo.Odor(i,1);
    end

    % Binarize TrialTrace and Flag Perturbations
    TrialTrace = Traces.Trial{i};
    TrialTrace(TrialTrace>0) = 1;

    if ~isempty(find(perturbedTrials==i))
        TrialTrace = TrialTrace*-1;
    end
    Traces.Trial{i} = TrialTrace;

end

%% Concatenate traces
whichTrials = 1:length(TrialInfo.TrialID);
traceOverlap = SampleRate*startoffset;
whichTraces{1} = 'Lever'; whichTraces{2} = 'Motor'; whichTraces{3} = 'Sniffs';
whichTraces{4} = 'Trial'; whichTraces{5} = 'Odor';
whichTraces{6} = 'Rewards'; whichTraces{7} = 'Timestamps';

for j = 1:size(whichTraces,2)
    temp = cellfun(@(x) ...
        x(1:end-traceOverlap), Traces.(whichTraces{j})(whichTrials), ...
        'UniformOutput', false);
    TracesOut.(whichTraces{j}) = {[cell2mat(temp(:)); ...
        Traces.(whichTraces{j}){whichTrials(end)}(end-traceOverlap+1:end,1)]};
end

% Sniffing specifc
% add a filtered sniff trace
TracesOut.SniffsFiltered{1}     = FilterThermistor(TracesOut.Sniffs{1});
% add a digital sniff trace
load(WhereSession,'CuratedSniffTimestamps');
if exist('CuratedSniffTimestamps','var')
    if size(CuratedSniffTimestamps,2) < 10
        CuratedSniffTimestamps(:,10) = 0;
    end
    LocationSniffs = TracesOut.SniffsFiltered{1}*nan;
    DigitalSniffs = TracesOut.SniffsFiltered{1}*0;
    for n = 1:size(CuratedSniffTimestamps)
        idx = CuratedSniffTimestamps(n,8:9);
        if CuratedSniffTimestamps(n,10) == -1
            DigitalSniffs(idx(1):idx(2)) = -1;
        else
            DigitalSniffs(idx(1):idx(2)) = 1;
        end
        if ~any(isnan(CuratedSniffTimestamps(:,4)))
            location = CuratedSniffTimestamps(n,4);
            LocationSniffs(idx(1):idx(2)) = location;
        end
    end

    TracesOut.SniffsDigitized{1} = DigitalSniffs;
    TracesOut.SniffsLocationed{1} = LocationSniffs;
end

% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% snity check for timestamp drops
if any(round(diff(Timestamps),3,'decimal')~=0.002)
    keyboard;
end

% calculate time adjustment
ephys = TTLs.Trial(1:numel(TrialInfo.Odor),2);
%behavior = TrialInfo.SessionTimestamps(:,2);
behavior = Timestamps(find((diff(abs(TracesOut.Trial{1})))==-1))';
myfit = fit(behavior,ephys-behavior,'poly1');

TimestampAdjust.ClosedLoop(2) = myfit.p2; 
TimestampAdjust.ClosedLoop(1) = myfit.p1; 

% if any(abs(TTLs.Trial(1:numel(TrialInfo.Odor),2) - (TrialInfo.SessionTimestamps(:,2) + TimestampAdjust.ClosedLoop))>0.04)
% TracesOut.Timestamps{1} = Timestamps + TimestampAdjust.ClosedLoop;
% sanity check for clock drift
if ~any(abs(TTLs.Trial(1:numel(TrialInfo.Odor),2) - (behavior + TimestampAdjust.ClosedLoop(2)))>0.04)
    % convert the behavior timestamps to OEPS base
    TracesOut.Timestamps{1} = (Timestamps + Timestamps*TimestampAdjust.ClosedLoop(1) + TimestampAdjust.ClosedLoop(2))';
else
    disp('clock drift in ephys and behavior files');
    keyboard;
end

%% tally odor valve timings with ephys TTLs
for x = 1:3
    odorvector = TracesOut.Odor{1};
    odorvector(odorvector~=x) = 0;
    odorTS = TracesOut.Timestamps{1}(find(diff(odorvector)));
    odorTS = reshape(odorTS,2,floor(numel(odorTS)/2))';
    
    if ~any(abs(odorTS(:,2)-TTLs.(['Odor',num2str(x)])(1:size(odorTS,1),2))>0.005)
        if any(abs(odorTS(:,1)-TTLs.(['Odor',num2str(x)])(1:size(odorTS,1),1))>0.005)
            f = find(abs(odorTS(:,1)-TTLs.(['Odor',num2str(x)])(1:size(odorTS,1),1))>0.005);
            for n = 1:numel(f)
                idx1 = find(TracesOut.Timestamps{1}==odorTS(f(n),1));
                idx2 = find(TracesOut.Timestamps{1}==odorTS(f(n),2));
                TracesOut.Odor{1}(idx1+1:idx2) = 0;
                [~,idx3] = min(abs(TracesOut.Timestamps{1}-TTLs.(['Odor',num2str(x)])(f(n),1)));
                TracesOut.Odor{1}(idx3:idx2) = x;
            end
        end
    else
       keyboard;
    end
end

%% for the manifold air and air(odor port)
AirVector = ~TracesOut.Odor{1}; % all periods when none of the odor ports are on
AirIdx = find(diff(AirVector));
AirTS  = TracesOut.Timestamps{1}(AirIdx);

if AirVector(1) == 1
    AirTS = vertcat(nan,AirTS);
end
if mod(numel(AirTS),2)
    AirTS = vertcat(AirTS,nan);
end
AirTS = reshape(AirTS,2,[])';
AirVector = (~AirVector)*1;

if ~any(abs(AirTS(2:end-1,1)-TTLs.Air(2:(size(AirTS,1)-1),1))>0.005) && ...
        ~any(abs(AirTS(2:end-1,2)-TTLs.Air(2:(size(AirTS,1)-1),2))>0.005)
    % check the first transition 
    if isnan(AirTS(1,1)) && isnan(TTLs.Air(1,1))
        if ~isequal(AirTS(1,2),TTLs.Air(1,2))
            [~,idx] = min(abs(TracesOut.Timestamps{1}-AirTS(1,2)));
            AirVector(1:idx) = -1; % assume all odors are off
            idx = find(TracesOut.Timestamps{1}<=TTLs.Air(1,2),1,'first');
            if ~isempty(idx)
                AirVector(1:idx+1) = 0;
            end
        end
    elseif ~isequal(AirTS(1,1),TTLs.Air(1,1))
        keyboard;
    end
    
    % check the last transition
    if abs(AirTS(end,1)-TTLs.Air(size(AirTS,1),1))<0.005
        [~,idx] = min(abs(TracesOut.Timestamps{1}-AirTS(end,1)));
        AirVector(idx+1:end) = 0; % assume Air stayed on
        idx = find(TracesOut.Timestamps{1}>=TTLs.Air(size(AirTS,1),2),1,'first');
        if ~isempty(idx)
            AirVector(idx:end) = -1; % assume Air stayed on
        end
    else
        keyboard;
    end
else
    keyboard;
end

manifoldVector = TracesOut.Odor{1};
manifoldVector(manifoldVector>0) = 1;
ManifoldIdx = find(diff(manifoldVector));
ManifoldIdx = reshape(ManifoldIdx,2,floor(numel(ManifoldIdx)/2))';
odorTS = TracesOut.Timestamps{1}(ManifoldIdx);
if ~any(abs(odorTS(:,1)-TTLs.AirManifold(1:size(odorTS,1),1))>0.005)
    for n = 1:size(ManifoldIdx,1)
        [~,idx] = min(abs(TracesOut.Timestamps{1}-TTLs.AirManifold(n,2)));
        manifoldVector(ManifoldIdx(n,1):idx) = 1;
    end
    TracesOut.Manifold{1} = manifoldVector;
else
    keyboard;
end

TracesOut.Odor{1} = AirVector;

%% Get the timestamp for Closed Loop End
SessionLength = ceil(TrialInfo.SessionTimestamps(end,2) + TimestampAdjust.ClosedLoop); % in seconds

%% for passive tuning part of the session
PassiveOut = MakePassiveSessionTraces(WhereSession);

%%
if savemode
    %% Organize data for Wolf separately
    path = fileparts(fileparts(fileparts(WhereSession)));
    if ~exist(fullfile(path,'WDW'),'dir')
        mkdir(fullfile(path,'WDW'));
    end
    [~,filename] = fileparts(WhereSession);
    filename = [filename,'.mat'];
    savepath = '/home/priyanka/Desktop';
    save(fullfile(savepath,'forWDW',filename),'TracesOut','PassiveOut','SingleUnits');
end

end
