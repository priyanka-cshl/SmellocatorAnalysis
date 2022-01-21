function [MyTrials, TrialSequence, ReplayTraces] = ParseTuningSession(FileName, PIDflag)
if nargin<2
    PIDflag = 0;
end

%% globals
global SampleRate; % = 500; % samples/second
global startoffset; % = 1; % seconds
global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips]

[MyData, MyParams, DataTags, TrialSequence] = LoadSessionData(FileName, 1, PIDflag); %#ok<ASGLU>

%% Get Trial ON-OFF timestamps
TrialColumn = MyData(:,6);
TrialColumn(TrialColumn~=0) = 1; % make logical
TrialOn = find(diff(TrialColumn)>0);
TrialOff =  find(diff(TrialColumn)<0)+1;

% account for cases where acquisition started/ended in between a trial
while TrialOn(1)>TrialOff(1)
    TrialOff(1,:) = [];
end
while TrialOn(end)>TrialOff(end)
    TrialOn(end,:) = [];
end

MyTrials = [NaN*ones(length(TrialOn),2) TrialOn TrialOff MyData(TrialOn,1) MyData(TrialOff,1)];
MyTrials(:,7) = MyTrials(:,6) - MyTrials(:,5); 
% MyTrials = [Location, OdorId, T-ON-idx, T-OFF-idx, T-ON-TS, T-OFF-TS, T-duration-TS]

%% Odor ON-OFF timestamps
OdorColumn = MyData(:,find(ismember(DataTags,'InRewardZone'))); % In reward zone column was used
OdorOn = find(diff(OdorColumn)>0);
OdorOff =  find(diff(OdorColumn)<0)+1;

while OdorOn(1)>OdorOff(1)
    OdorOff(1,:) = [];
end
while OdorOn(end)>OdorOff(end)
    OdorOn(end,:) = [];
end
OdorValve = [OdorOn OdorOff];

MotorColumn = MyData(:,find(ismember(DataTags,'Motor')));

% Fill up the Trials Table with odor on-off timestamps
for i = 1:size(MyTrials,1)
    % add 4 more columns to MyTrials
    % col 8,9 - odor valve ON and OFF idx, col 10,11 are the timestamps
    if ~isempty(intersect(find(OdorValve(:,1)>MyTrials(i,3)),find(OdorValve(:,1)<MyTrials(i,4))))
        MyTrials(i,8:9) = OdorValve(intersect(find(OdorValve(:,1)>MyTrials(i,3)),find(OdorValve(:,1)<MyTrials(i,4))),:);
        MyTrials(i,10:11) = MyData(MyTrials(i,8:9),1); % timestamp
    else
        MyTrials(i,8:11) = 0;
    end
    % get the motor position
    x2 = MyTrials(i,4);
    x1 = max(MyTrials(i,3),MyTrials(i,8));
    MyTrials(i,1) = round(mode(MotorColumn(x1:x2,1)));
end

TrialSequence = vertcat([NaN NaN], TrialSequence(1:end-1,:));

% if there are any passive replay trials - extract the traces
if any(TrialSequence(:,1)==999)
    
    % which trials
    x = find(MyTrials(:,7)>sum(MyParams(3:8))/1000); % these trials will be longer
    
    % get column IDs
    LeverCol = find(cellfun(@isempty,regexp(DataTags,'Lever'))==0);
    MotorCol = find(cellfun(@isempty,regexp(DataTags,'Motor'))==0);
    if ~isempty(find(cellfun(@isempty,regexp(DataTags,'thermistor'))==0))
        RespCol = find(cellfun(@isempty,regexp(DataTags,'thermistor'))==0);
    else
        RespCol = find(cellfun(@isempty,regexp(DataTags,'respiration'))==0);
    end
    LickCol = find(cellfun(@isempty,regexp(DataTags,'Licks'))==0);
    TrialCol = find(cellfun(@isempty,regexp(DataTags,'TrialON'))==0);
    RewardCol = find(cellfun(@isempty,regexp(DataTags,'Rewards'))==0);
    
    % extract traces
    for i = 1:numel(x)
        thisTrial = x(i); % trial #
        start_idx = MyTrials(thisTrial,3) - startoffset*SampleRate; % 1 sec preceding trial start
        stop_idx = MyTrials(thisTrial+1,3) -1; % upto next trial start
        
        ReplayTraces.Lever(i)     = { MyData(start_idx:stop_idx, LeverCol) };
        ReplayTraces.Motor(i)     = { MyData(start_idx:stop_idx, MotorCol) };
        ReplayTraces.Sniffs(i)    = { MyData(start_idx:stop_idx, RespCol) };
        ReplayTraces.Licks(i)     = { MyData(start_idx:stop_idx, LickCol) };
        ReplayTraces.Trial(i)     = { MyData(start_idx:stop_idx, TrialCol) };
        ReplayTraces.Rewards(i)   = { MyData(start_idx:stop_idx, RewardCol) };
        
        % extract the timestamps incase it can't be reconstructed from indices
        ReplayTraces.Timestamps(i)  = { MyData(start_idx:stop_idx, 1) };
    end
else
    ReplayTraces = [];
end

end