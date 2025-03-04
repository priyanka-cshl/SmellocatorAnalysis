function [TrialAligned, TrialInfo, ReplayAligned, ReplayInfo, TuningAligned, TuningInfo, AllUnits] = PreprocessSpikesAndSniffs(MySession)
% general script to take the output of PreprocessSmellocatorData
% and organize the ephys data by aligning to trials
% and also to trial-aligned sniffs

% initialize some globals - they are used in replay processing
global SampleRate; % = 500; % samples/second
global startoffset; % = 1; % seconds

%% 1: Load the preprocessed session
%   Load the relevant behavior and ephys related variables
load(MySession, 'TrialInfo', 'TTLs', 'SampleRate', 'startoffset');
%   Keep track of SessionPath - useful when passing around to other functions
TrialInfo.SessionPath = MySession;

%% 2: Get Relevant Trial Aligned Events for all trials
TrialAligned.EventLegends = ...
    {'OdorStart', 'Reward', 'TrialOFF', 'First TZ Entry', 'PerturbationStart'};

nTrials         = size(TrialInfo.TrialID,2);

% Somethings have already been computed during Behavior+Ephys Preprocessing
OdorStart       = TTLs.Trial(1:nTrials,4); % already computed during Ephys Preprocessing in PreprocessSmellocator
TrialOFF        = TTLs.Trial(1:nTrials,3); % already computed during Ephys Preprocessing in PreprocessSmellocator

% Some need to be computed fresh
% Perturbation start : time-point equivalent to perturbartion start in control trials
[whichPerturbations, ia, ic] = ...
    unique(TrialInfo.Perturbation(find(~cellfun(@isempty, TrialInfo.Perturbation(:,1))),1));
perturbation_counts = accumarray(ic,1);
% ignore any perturbation that occurs just once
whichPerturbations(find(perturbation_counts<2),:) = [];
% find only the within trial perturbations
if numel(find(~strcmp(whichPerturbations,'OL-Template'))) == 1
    PerturbationTag = whichPerturbations(find(~strcmp(whichPerturbations,'OL-Template')));
    PerturbationTag = PerturbationTag{1};
else
    disp('Multiple different Perturbations found!');
    keyboard;
end

[PerturbationEvents] = FindPutativePerturbationStart(TrialInfo, PerturbationTag);
% for halts:    [nx3]:[HaltStart HaltEnd HaltLocation]
% for offsets:  [nx3]:[[OffsetStart FeedbackResume OffsetLocation]
% for reverals: [nx1]:[[1 if reversed, 0 otherwise]
% n = trials

if ~isempty(PerturbationEvents)
    load(TrialInfo.SessionPath,'SampleRate');
    PerturbStart = PerturbationEvents(:,1)/SampleRate;
else
    PerturbStart    = NaN*(ones(nTrials,1));
end

FirstTZEntry    = NaN*(ones(nTrials,1));
Reward          = NaN*(ones(nTrials,1));

for whichtrial = 1:nTrials % every trial
    % First TZ entry
    if ~isempty(TrialInfo.InZone{whichtrial}) % Periods of 'being within' Target
        thisTrialEntries = TrialInfo.InZone{whichtrial};
        thisTrialEntries(:,3) = thisTrialEntries(:,2)-thisTrialEntries(:,1);
        if ~isempty(thisTrialEntries(find(thisTrialEntries(:,3)>=0.1,1,'first'),1))
            FirstTZEntry(whichtrial,1) = thisTrialEntries(find(thisTrialEntries(:,3)>=0.1,1,'first'),1);
        else
            FirstTZEntry(whichtrial,1) = thisTrialEntries(1,1);
        end
    end
    
    % Reward: overkill because typically its the end of the trial
    x1 = TTLs.Trial(whichtrial,1); % trial start
    x2 = TTLs.Trial(whichtrial,2) + 0.2; % next trial cannot start sooner than that
    rewardTS = TTLs.Reward(find((TTLs.Reward(:,1)>x1)&(TTLs.Reward(:,1)<x2)),1);
    if ~isempty(rewardTS)
        Reward(whichtrial,1) = rewardTS(1) - x1;
    end
end

TrialAligned.Events = [OdorStart Reward TrialOFF FirstTZEntry PerturbStart];

%% 3: Minor additions to TrialInfo
%   Keep track of SessionPath - useful when passing around to other functions
% TrialInfo.SessionPath = MySession;

% transpose some fields - just to keep stuff consistent (rows = trials)
TrialInfo.TrialID   = TrialInfo.TrialID';
TrialInfo.Offset    = TrialInfo.Offset';
TrialInfo.Reward    = TrialInfo.Reward';

%   Add a field (FirstEntry) to mark time-point after which the lever is largely in the Target Zone
TrialInfo.TargetEntry = FirstTZEntry;

% If there are Rule reversals then relabel control trials after a flip block
if any(find(strcmp(TrialInfo.Perturbation,'RuleReversal')))
    count = 0;
    y = find(strcmp(TrialInfo.Perturbation(:,1),'RuleReversal'));
    if numel(y)==1
        TrialInfo.Perturbation(y,1) = {[]};
    else
        if any(diff(y)~=1)
            x1 = y(find(diff(y)~=1))+1;
            x2 = y(find(diff(y)~=1)+1)-1;
            count = count + 1;
            TrialInfo.Perturbation(x1:x2,1) = {['RuleNormal',num2str(count)]};
            x3 = x2 + 1;
            TrialInfo.Perturbation(x3:y(end),1) = {['RuleReversal',num2str(count)]};
        end
        if TrialInfo.TrialID(end)>y(end)
            count = count + 1;
            TrialInfo.Perturbation((y(end)+1):end,1) = {['RuleNormal',num2str(count)]};
        end
    end
end

%% 4: Process replays
%   Parse Replays into subtrials and then contsruct equivalents of
%   TrialInfo (ReplayInfo) and TrialAlignedEvents (ReplayAlignedEvents)
load(TrialInfo.SessionPath, 'Traces', 'PassiveReplayTraces', 'ReplayTTLs');
if any(find(strcmp(TrialInfo.Perturbation(:,1),'OL-Replay'))) || ...
        ~isempty(PassiveReplayTraces)
    
    % process the replays
    OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);
    
    % parse replays into subtrials
    ReplayInfo = ParseReplaysToSubtrials(OpenLoop,TrialInfo,ReplayTTLs,TTLs);
    
    ReplayAligned.EventsLegends = TrialAligned.EventLegends;
    ReplayAligned.Events = TrialAligned.Events(ReplayInfo.TemplateInfo(:,2),:);
else
    ReplayInfo = [];
    ReplayAligned.Events = [];
    ReplayAligned.EventsLegends = [];
end

%% 5: Process Tuning
load(TrialInfo.SessionPath, 'Tuning*');
% loads TuningTTLs and Tuningextras;
if ~exist('Tuningextras')
    TuningSettings = [];
    TuningSequence = [];
else
    TuningSettings = Tuningextras.sessionsettings;
    TuningSequence = Tuningextras.sequence;
end

TuningInfo = ParseTuningSubtrials(TuningTTLs, TuningSequence, TuningSettings);
TuningAligned.Events = [];
TuningAligned.EventsLegends = [];

%% 6: Process Sniffs
% for the closed loop part
load(TrialInfo.SessionPath,'CuratedSniffTimestamps');
if exist('CuratedSniffTimestamps','var')
    SniffTS = CuratedSniffTimestamps(:,1:4);
else
    load(TrialInfo.SessionPath,'SniffTS');
end
% for the passive part
load(TrialInfo.SessionPath,'CuratedPassiveSniffTimestamps');
if exist('CuratedPassiveSniffTimestamps','var')
    SniffTS_passive = CuratedPassiveSniffTimestamps(:,1:4);
else
    load(TrialInfo.SessionPath,'SniffTS_passive');
end

% Sniffs are in behavior timebase
load(TrialInfo.SessionPath,'TimestampAdjust');
if isfield(TimestampAdjust,'ClosedLoop')
    if length(TimestampAdjust.ClosedLoop)==1
        if any(abs(TTLs.Trial(1:numel(TrialInfo.Odor),2) - (TrialInfo.SessionTimestamps(:,2) + TimestampAdjust.ClosedLoop))>0.04)
            disp('clock drift in ephys and behavior files');
            keyboard;
        end
        AllSniffs = [ (SniffTS(:,1:3) + TimestampAdjust.ClosedLoop)  SniffTS(:,4); ...
        (SniffTS_passive(:,1:3) + TimestampAdjust.Passive) SniffTS_passive(:,4) ];
    end
else
    % calculate the adequate time correction
    ephys = TTLs.Trial(1:numel(TrialInfo.Odor),2);
    behavior = TrialInfo.SessionTimestamps(:,2);
    myfit = fit(behavior,ephys-behavior,'poly1');

    TimestampAdjust.ClosedLoop(2) = myfit.p2; 
    TimestampAdjust.ClosedLoop(1) = myfit.p1; 
    
    AllSniffs = [ (SniffTS(:,1:3) + SniffTS(:,1:3).*TimestampAdjust.ClosedLoop(1) + TimestampAdjust.ClosedLoop(2))  SniffTS(:,4); ...
        (SniffTS_passive(:,1:3) + TimestampAdjust.Passive) SniffTS_passive(:,4) ];
end

% For Sniff Alingment to Trials
% need a 4 col TrialTimes matrix [nTrials x 4]
% [TrialStart TrialEnd PrevTrialOFF NextTrialStart];
% cols 3 and 4 basically create the window bounds around trialON (col1)

BufferWindow = 1; % in seconds,

% Get TrialTimes matrix for the close-loop session (in OEPS base)
TrialTimes = TTLs.Trial(1:nTrials,1:2); % TrialOn and TrialOFF in OEPS base
if size(TrialTimes,2) == 2
    % this will happen when passing trialtimes from the closed loop session
    % compute the previous trial off and next trial ON times
    TrialTimes(1,3)         = 0; % session start
    TrialTimes(2:end,3)     = TrialTimes(1:end-1,2);
    TrialTimes(1:end-1,4)   = TrialTimes(2:end,1);
    TrialTimes(end,4)       = TrialTimes(end,2) + BufferWindow; % 1 sec after the last trial
end

% Align sniffs by trials
[TrialAligned.Sniffs, TrialAligned.RefTS] = GetWindowAlignedSniffs(TrialTimes, AllSniffs, TrialInfo); % both in OEPS time base
% type help GetWindowAlignedSniffs to check notes on the output
TrialAligned.SniffLegends = ...
    {'Trial/Event ID', 'Odor/Air State', 'EventRefSniffIndex', 'SniffType', ...
     'InhStart_Previous', 'InhEnd_Previous', 'InhStart_Current', 'InhEnd_Current', 'InhStart_Next', 'InhEnd_Next', 'InhStart_NextNext', ...
     'MotorLocation_Previous', 'MotorLocation_Current', 'MotorLocation_Next'};

% Align sniffs by replay subtrials
% includes active/passive replays and passive perturbation replays
% replay trial times in ReplayInfo and are already in OEPS base
if ~isempty(ReplayInfo)
    [ReplayAligned.Sniffs, ReplayAligned.RefTS] = GetWindowAlignedSniffs(ReplayInfo.trialtimes,AllSniffs,ReplayInfo);
    ReplayAligned.SniffLegends = TrialAligned.SniffLegends;
end

% Align sniffs by tuning subtrials
[TuningAligned.Sniffs, TuningAligned.RefTS] = GetWindowAlignedSniffs(TuningInfo.trialtimes,AllSniffs,TuningInfo);
TuningAligned.SniffLegends = TrialAligned.SniffLegends;
TuningAligned.SniffLegends{4} = 'Ignore';

%% Process Spikes
% simplify stuff here
% just get the raw timestamps and tetrode + channel info for each unit
load(TrialInfo.SessionPath,'SingleUnits');
for whichunit = 1:size(SingleUnits,2)
    AllUnits.Spikes{whichunit}      = SingleUnits(whichunit).spikes; % raw timestamps in OEPS base
    AllUnits.ChannelInfo(whichunit,1:2) = [SingleUnits(whichunit).tetrode SingleUnits(whichunit).id]; % tetrode and phy cluster ID
    AllUnits.SpikeAmps{whichunit} = SingleUnits(1).spikescaling(find(SingleUnits(1).clusterscalingorder == SingleUnits(whichunit).id));
end

end