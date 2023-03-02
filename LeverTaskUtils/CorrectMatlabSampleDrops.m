% organize the session data into a cell array of trials
function [Trial] = ...
    CorrectMatlabSampleDrops(MyData, MySettings, DataTags)

% Check if there were any sample drop issues between digital and analog data
% this issue was observed for batch K on rig2

% Outputs: 
% Trial.Indices = trial start, trial stop and duration (in indices)
% Trial.TimeStamps = trial start, trial stop and duration (in timestamps)
% Trial.Offsets = trial start offset and odor start offset, from digital trial start (in indices)

%% globals
global SampleRate; % = 500; % samples/second
global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips-home sensor, motor skips-TZ mismatch]

%% Get Column IDs
TrialCol = find(cellfun(@isempty,regexp(DataTags,'TrialON'))==0);
LeverCol = find(cellfun(@isempty,regexp(DataTags,'Lever'))==0);

%% Get Trial ON-OFF timestamps
TrialColumn = MyData(:,TrialCol);
TrialColumn(TrialColumn~=0) = 1; % make logical
TrialOn = find(diff([0; TrialColumn])>0); % this should work better in case samples were dropped in NI : PG 23/03/02
% TrialOn = find(diff([0; TrialColumn])>0) -1;
TrialOff =  find(diff(TrialColumn)<0)+1;

% account for cases where acquisition started/ended in between a trial
while TrialOn(1)>TrialOff(1)
    TrialOff(1,:) = [];
end
while TrialOn(end)>TrialOff(end)
    TrialOn(end,:) = [];
end 

Trial.Indices = [TrialOn TrialOff (TrialOff-TrialOn)];

% if acq started when trial was already high - TrialOn(1) == 0 : PG 23/03/02 - no longer
% valid - TrialOn(1) cannot be less than 1 - next line still works though
Trial.TimeStamps(find(TrialOn),1) = MyData(TrialOn(find(TrialOn)),1); 

Trial.TimeStamps(:,2) = MyData(TrialOff,1);
Trial.TimeStamps(:,3) = Trial.TimeStamps(:,2) - Trial.TimeStamps(:,1);

if Trial.TimeStamps(1,1) == 0 % PG 23/03/02 - no longer valid - TrialOn(1) cannot be less than 1
    Trial.Indices(1,1) =  1;
end

% check for timestamp drops
TimeStamp_drops = find(abs(diff(MyData(:,1)))>((1/SampleRate)+0.0001)) + 1;
if ~isempty(TimeStamp_drops)
    errorflags(2) = 1;
    % these timestamp drops seem to happen only in chunks with Trial Starts
    % right at the end - check this
    if isequal(TimeStamp_drops,intersect(TrialOn,TimeStamp_drops))
        disp(['Warning: timestamps dropped at ', num2str(numel(TimeStamp_drops)), ' trial starts']);
    else
        disp('Warning: timestamps dropped - not just at trial starts - ignoring those');
        keyboard;
    end
    for i = 1:numel(TimeStamp_drops)
        f = find(TrialOn==TimeStamp_drops(i));
        if ~isempty(f)
            Trial.TimeStamps(f,1) = Trial.TimeStamps(f,1) - (1/SampleRate);
            Trial.TimeStampDrops(f,1) = 1;
        end
    end
    % recompute trial durations
    Trial.TimeStamps(:,3) = Trial.TimeStamps(:,2) - Trial.TimeStamps(:,1);
else
   disp('No timestamps dropped');
   Trial.TimeStampDrops(1:numel(TrialOn),1) = 0;
end

% for trial On detection 
LeverThresh = median(MySettings(:,11));
trialflag = [];

%% Process lever snippets trial-by-trial
for thisTrial = 1:size(Trial.Indices,1)
    
    if TrialOn(thisTrial)
    
        % extract continuous traces for lever
        % extract from previous trial OFF and upto this trial start
        if thisTrial == 1 % all except first trial
            LastTrialIdx = 1;
        else
            LastTrialIdx = TrialOff(thisTrial-1);
        end
        
        thisTrialIdx = TrialOn(thisTrial);
        
        %% Timestamps for Trial ON - reconstructed from Lever trace
        % Odor ON timestamp
        LeverSnippet = MyData(LastTrialIdx:thisTrialIdx, LeverCol);
        LeverTemp = LeverSnippet;
        % assume a fix threshold of 4.8
        LeverTemp(LeverSnippet<=LeverThresh) = 0;
        LeverTemp(LeverTemp>0) = 1;
        % get all initiation hold periods
        Initiations = [find(diff([0; LeverTemp; 0])==1) find(diff([0; LeverTemp; 0])==-1)-1];
        TS = MyData(LastTrialIdx:thisTrialIdx, 1);
        if Trial.TimeStampDrops(thisTrial) && Initiations(end) == (numel(LeverSnippet)-1)
            Initiations(end) = Initiations(end)+1;
        end
        Initiations_durations = 1000*diff([TS(Initiations(:,1),1) TS(Initiations(:,2),1)],1,2); % in ms

        if ~isempty(Initiations)
            % Odor ON timestamp - find the first initiation period > trigger hold
            TriggerHold = MySettings(...
                                find(MySettings(:,1)<=Trial.TimeStamps(thisTrial,1),1,'last') ...
                                    ,13); % in msec
            OdorStart = find(Initiations_durations>=TriggerHold, 1, 'last');
            if isempty(OdorStart)
                if size(Initiations,1) > 1
                    % no initiation bout longer than TriggerHold (bug)
                    % take the largest bout
                    % [~,OdorStart] = max(diff(Initiations,1,2));
                    % changed for batch Q: I don't know why we would ever pick anything but the last one 
                    OdorStart = size(Initiations,1);
                else
                    % take the initiation bout thats available
                    OdorStart = 1;
                end
                %disp(['Trial ',num2str(thisTrial),': No valid Initiations found!']);
                trialflag(thisTrial) = -1;
            elseif OdorStart < size(Initiations,1) % changed for batch Q:
                % I don't know why we would ever pick anything but the last one 
                OdorStart = size(Initiations,1);
                disp(thisTrial);
                trialflag(thisTrial) = -1;
            end
            
            % Trial ON timestamp - find the first initiation period > trigger hold
            TrialStartOffsets(thisTrial,1) = Initiations(OdorStart,2) - numel(LeverSnippet);
            %OdorStartOffsets(thisTrial,1) = Initiations(OdorStart,1) + TriggerHold - numel(LeverSnippet);
            OdorStartOffsets(thisTrial,1) = TS(Initiations(OdorStart,1)) - TS(end) + TriggerHold/1000;
        else
            %TrialStartOffsets(thisTrial,1) = NaN;
            OdorStartOffsets(thisTrial,1) = NaN;
        end
    else
        %TrialStartOffsets(thisTrial,1) = NaN;
        OdorStartOffsets(thisTrial,1) = NaN;
    end
    
end

% Some offsets are miscalculated - as being zero due to the threshold not
% being exactly == 4.8V - this HACK is an attempt to eliminate those
% do this only if there are any sample drops in the first place
if numel(find(TrialStartOffsets<-5))>=5 && isempty(find(cellfun(@isempty,regexp(DataTags,'thermistor'))==0))
    spurious_offsets = 1;
    while spurious_offsets
        % pick the last offset that's non-zero
        f = find(TrialStartOffsets==0,1,'last');
        % check if the offsets in the trial before and after are
        % significantly different from 0
        if f == numel(TrialStartOffsets) && abs(TrialStartOffsets(f-1))>5
            % force it to be the same as the trial before
            TrialStartOffsets(f) = TrialStartOffsets(f-1);
            OdorStartOffsets(f) = NaN;
            trialflag(f) = -2;
        elseif f == 1
            % makes no sense to have sample drops on the first trial itself
            % ignore this offset
            TrialStartOffsets(f) = 0;
            OdorStartOffsets(f) = 0;
        elseif abs(TrialStartOffsets(f-1))>5 && abs(TrialStartOffsets(f+1))
            % force it to be the same as the trial before
            TrialStartOffsets(f) = TrialStartOffsets(f-1);
            OdorStartOffsets(f) = NaN;
            trialflag(f) = -2;
        else
            spurious_offsets = 0;
            break;
        end
    end
end

% ignore any very large offsets
TrialStartOffsets(TrialStartOffsets<-200) = NaN;

% check if there's just one or two offsets - Rig I
% in that case, its just due to noise in the lever signal
% ignore all offsets
foo = TrialStartOffsets;
foo(isnan(foo)) = 0;
% threshold
foo(foo>-5) = 0;
foo(foo<0) = 1;
% any contiguous stretch of 5 offsets?
x = [find(diff([0; foo; 0])==1) find(diff([0; foo; 0])==-1)-1];
Trial.Offsets = [any(abs(diff(x,1,2))>5)*TrialStartOffsets OdorStartOffsets];

if ~any(abs(diff(x,1,2))>5)
    disp('No Samples dropped between digital and analog');
else
    disp('Warning: Samples dropped between digital and analog');
    errorflags(1) = 1;
end

end