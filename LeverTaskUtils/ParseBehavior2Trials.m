% organize the session data into a cell array of trials
function [Traces, TrialInfo] = ...
    ParseBehavior2Trials(MyData, MySettings, DataTags, Trial)

% Extract traces starting from 1 sec before trial start to next trial start
% First trial - 1 sec before trial start (or from session start if smaller)
% Last trial - 1 sec after this trial stop (or session stop if smaller)
% traces will overlap by 1 sec

%% globals
global SampleRate; % = 500; % samples/second
global startoffset; % = 1; % seconds
global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips]
global TargetZones;

%% get list of target zones used
TargetZones = unique(MySettings(:,18:20),'rows');
% bug handler - ignore spurious target zones
if size(TargetZones,1)>12
    disp('WARNING: session contains more than 12 target zones');
    TargetZonesTemp = TargetZones;
    f = find(round(TargetZones(:,1) - TargetZones(:,3))~=1);
    TargetZones(f,:) = [];
end

%% get list of fake target zones used
FakeTargetZones = unique(MySettings(:,26:28),'rows');
% bug handler - ignore spurious fake target zones
foo = FakeTargetZones;
foo(:,2) = foo(:,2) - foo(:,1);
if ~isempty(find((foo(:,2)==0)&(foo(:,1)<20)&(foo(:,1)>0)))
    disp('WARNING: session contains buggy fake zones');
    FakeTargetZones(find((foo(:,2)==0)&(foo(:,1)<20)&(foo(:,1)>0)),:) = [];
end

%% Get Column IDs
TrialCol = find(cellfun(@isempty,regexp(DataTags,'TrialON'))==0);
LeverCol = find(cellfun(@isempty,regexp(DataTags,'Lever'))==0);
MotorCol = find(cellfun(@isempty,regexp(DataTags,'Motor'))==0);
LickCol = find(cellfun(@isempty,regexp(DataTags,'Licks'))==0);
RewardCol = find(cellfun(@isempty,regexp(DataTags,'Rewards'))==0);
TZoneCol = find(cellfun(@isempty,regexp(DataTags,'InTargetZone'))==0);
RZoneCol = find(cellfun(@isempty,regexp(DataTags,'InRewardZone'))==0);
if ~isempty(find(cellfun(@isempty,regexp(DataTags,'thermistor'))==0))
    RespCol = find(cellfun(@isempty,regexp(DataTags,'thermistor'))==0);
else
    RespCol = find(cellfun(@isempty,regexp(DataTags,'respiration'))==0);
end
PerturbationCol(1) = find(cellfun(@isempty,regexp(DataTags,'WhichPerturbation'))==0);
PerturbationCol(2) = find(cellfun(@isempty,regexp(DataTags,'PerturbationValue'))==0);
%% Get Trial ON-OFF indices and timestamps
TrialOn  = Trial.Indices(:,1);
TrialOff = Trial.Indices(:,2);
TrialOffsets = Trial.Offsets(:,1); % correct for the trial offset (sample drops on digital channel)
OdorOffsets = Trial.Offsets(:,2);

%% Crunch data trial-by-trial
for thisTrial = 1:numel(TrialOn)
    % store original trial ID - some trials may get deleted later because of weird target zones
    TrialInfo.TrialID(thisTrial) = thisTrial;
    thisTrialOffset = TrialOffsets(thisTrial); % this will be zero if there were no digital-analog sample drops
    
    if ~isnan(thisTrialOffset)
        % extract continuous traces for lever, motor position, licks and sniffs
        
        % correction factor
        TrialInfo.Offset(thisTrial) = thisTrialOffset;
        
        start_idx = TrialOn(thisTrial) - startoffset*SampleRate; % 1 sec preceding trial start
        start_idxCorrected = start_idx + thisTrialOffset;
        if thisTrial == 1 % exception handler for the first trial
            start_idx = max(1,start_idx);
            start_idxCorrected = max(1,start_idxCorrected);
            LastTrialIdx = start_idx;
        end
        
        if thisTrial < length(TrialOn)
            stop_idx = TrialOn(thisTrial+1) -1;
            stop_idxCorrected = stop_idx + thisTrialOffset;
        else  % exception handler for the last trial
            stop_idx = TrialOff(thisTrial) + startoffset*SampleRate;
            stop_idxCorrected = stop_idx + thisTrialOffset;
            stop_idx = min(stop_idx,size(MyData,1));
            stop_idxCorrected = min(stop_idxCorrected,size(MyData,1));
        end                
        
        %% Extract traces
        % Analog - use corrected indices - analog channnel dropped samples
        Traces.Lever(thisTrial) = { MyData(start_idxCorrected:stop_idxCorrected, LeverCol) };
        Traces.Motor(thisTrial) = { MyData(start_idxCorrected:stop_idxCorrected, MotorCol) };
        Traces.Sniffs(thisTrial) = { MyData(start_idxCorrected:stop_idxCorrected, RespCol) };
        
        Traces.Licks(thisTrial) = { MyData(start_idx:stop_idx, LickCol) };
        Traces.Trial(thisTrial) = { MyData(start_idx:stop_idx, TrialCol) };
        Traces.Rewards(thisTrial) = { MyData(start_idx:stop_idx, RewardCol) };
        
        % start and stop indices of the extracted trace - w.r.t to the session
        TrialInfo.TraceIndices(thisTrial,:) = [start_idx stop_idx start_idxCorrected stop_idxCorrected];
        TrialInfo.TraceDuration(thisTrial,1) = (diff([start_idx stop_idx]) + 1)/SampleRate;
        Traces.Timestamps(thisTrial) = { MyData(start_idx:stop_idx, 1) };
        % extract the timestamps if it can't be reconstructed from indices
        if errorflags(2) % timestamps were dropped
            Traces.TimestampsAnalog(thisTrial) = { MyData(start_idxCorrected:stop_idxCorrected, 1) };
        end
        
        % handing timestamp drops
        TrialInfo.TimeStampsDropped(thisTrial) = Trial.TimeStampDrops(thisTrial);
        if TrialInfo.TimeStampsDropped(thisTrial)
            % find the location where missing samples need to be patched in
           f = find(abs(diff(Traces.Timestamps{thisTrial}))>SampleRate^-1+0.0001);
           samples_missing = int16((Traces.Timestamps{thisTrial}(f+1) - Traces.Timestamps{thisTrial}(f))*SampleRate - 1);
           timestamps_missing = linspace(Traces.Timestamps{thisTrial}(f),Traces.Timestamps{thisTrial}(f+1),samples_missing+2);
           timestamps_missing(:,1) = [];
           timestamps_missing(:,end) = [];
           % pad NaNs into the Traces
           Traces.Lever{thisTrial} = [Traces.Lever{thisTrial}(1:f); Inf(samples_missing,1); Traces.Lever{thisTrial}(f+1:end)];
           Traces.Motor{thisTrial} = [Traces.Motor{thisTrial}(1:f); Inf(samples_missing,1); Traces.Motor{thisTrial}(f+1:end)];
           Traces.Sniffs{thisTrial} = [Traces.Sniffs{thisTrial}(1:f); Inf(samples_missing,1); Traces.Sniffs{thisTrial}(f+1:end)];
           Traces.Licks{thisTrial} = [Traces.Licks{thisTrial}(1:f); Inf(samples_missing,1); Traces.Licks{thisTrial}(f+1:end)];
           Traces.Rewards{thisTrial} = [Traces.Rewards{thisTrial}(1:f); Inf(samples_missing,1); Traces.Rewards{thisTrial}(f+1:end)];
           % making an assumption here
           Traces.Trial{thisTrial} = [Traces.Trial{thisTrial}(1:f); zeros(samples_missing-1,1); 1; Traces.Trial{thisTrial}(f+1:end)];
           Traces.Timestamps{thisTrial} = [Traces.Timestamps{thisTrial}(1:f); timestamps_missing'; Traces.Timestamps{thisTrial}(f+1:end)];
           Traces.TimestampsAnalog{thisTrial} = [Traces.TimestampsAnalog{thisTrial}(1:f); timestamps_missing'; Traces.TimestampsAnalog{thisTrial}(f+1:end)];
        end
        
        %% Extract Trial Timestamps
        % w.r.t SessionStart (to go back to raw data if needed)
        thisTrialIdx = [TrialOn(thisTrial) TrialOff(thisTrial)]; %uncorrected
        thisTrialIdxCorrected = thisTrialIdx + thisTrialOffset;        
        TrialInfo.SessionIndices(thisTrial,:) = [thisTrialIdx thisTrialIdxCorrected]; % actual session indices - to go back to raw data
        TrialInfo.SessionTimestamps(thisTrial,:) = MyData([thisTrialIdx thisTrialIdxCorrected],1); % actual timestamps of trial start and end
        % w.r.t TrialStart
        TrialInfo.TimeIndices(thisTrial,:) = max(1,thisTrialIdx - start_idx); % when session starts with trial ON - this value becomes zero
        TrialInfo.Timestamps(thisTrial,:) = MyData(thisTrialIdx,1) - MyData(start_idx,1); % in seconds
        %TrialInfo.Duration(thisTrial,1) = (diff(thisTrialIdx) + 1)/SampleRate; % in seconds
        TrialInfo.Duration(thisTrial,1) = diff(TrialInfo.Timestamps(thisTrial,:));
        
        %% Which odor
        TrialInfo.Odor(thisTrial,1) = mode(MyData(thisTrialIdx(1):thisTrialIdx(2),TrialCol));
        % Odor ON timestamp (from the InRewardZone column - enocdes Odor ON before trialstart - see GUI)
        thisTrialInZone = find(diff(MyData(LastTrialIdx:thisTrialIdx(1), RZoneCol))==-1);
        if TrialInfo.TimeStampsDropped(thisTrial) && thisTrialInZone(end) == (numel(LastTrialIdx:thisTrialIdx(1)) - 1)
                thisTrialInZone = []; % cannot determine accurately when the Odor went ON - transition samples are missing
        end

        if ~isempty(thisTrialInZone)
            %TrialInfo.OdorStart(thisTrial,1) = LastTrialIdx + thisTrialInZone(end) - thisTrialIdx(1); % odor start idx w.r.t trial start
            TrialInfo.OdorStart(thisTrial,1) = MyData(LastTrialIdx + thisTrialInZone(end) - 1,1) - MyData(thisTrialIdx(1),1); % odor start idx w.r.t trial start
        else
            TrialInfo.OdorStart(thisTrial,1) = OdorOffsets(thisTrial);
        end
        
        %% Timestamps for Odor ON and Trial ON - reconstructed from Lever trace
        TrialInfo.OdorStart(thisTrial,2) = OdorOffsets(thisTrial);
        %TrialInfo.OdorStart(thisTrial,:) = TrialInfo.OdorStart(thisTrial,:)/SampleRate;
        
        %% Which TargetZone
        if ~isempty(find(TargetZones(:,1) == mode(MyData(thisTrialIdx(1):thisTrialIdx(2),2)),1))
            TrialInfo.TargetZoneType(thisTrial,1) = ...
                find(TargetZones(:,1) == mode(MyData(thisTrialIdx(1):thisTrialIdx(2),2)),1);
        else
            % bug handler - for spurious target zones
            thiszonetarget = TargetZonesTemp(find(TargetZonesTemp(:,1) == mode(MyData(thisTrialIdx(1):thisTrialIdx(2),2)),1),2);
            TrialInfo.TargetZoneType(thisTrial,1) = find(TargetZones(:,2) == thiszonetarget);
        end
        
        %% TF : odor starts from left or right?
        % check the motor position at trialstart - 10 samples before trial start
        % to verify if the transfer function was inverted in this trial
        if thisTrial>1
            TrialInfo.TransferFunctionLeft(thisTrial,1) = (MyData(TrialOn(thisTrial)-1, MotorCol)>0);
        else
            TrialInfo.TransferFunctionLeft(thisTrial,1) = 0;
        end
        
        %% Reward timestamps
        thisTrialRewards = find(diff(MyData(start_idx:stop_idx,RewardCol))==1); % indices w.r.t. to trace start
        thisTrialRewards = thisTrialRewards/SampleRate; % convert to seconds
        % force the reward time stamps that were before trial start to be -ve
        thisTrialRewards(thisTrialRewards < TrialInfo.Timestamps(thisTrial,1)) = ...
            -1*thisTrialRewards(thisTrialRewards < TrialInfo.Timestamps(thisTrial,1));
        if ~isempty(thisTrialRewards)
            TrialInfo.Reward(thisTrial) = { thisTrialRewards };
            TrialInfo.Success(thisTrial,1) = any(thisTrialRewards>0); % successes and failures
        else
            TrialInfo.Reward(thisTrial) = { [] };
            TrialInfo.Success(thisTrial,1) = 0; % successes and failures
        end
        
        %% Calculate all stay times (in the target zone)
        thisTrialInZone = [find(diff([0;MyData(TrialOn(thisTrial):TrialOff(thisTrial), TZoneCol)])==1) ...
            find(diff([MyData(TrialOn(thisTrial):TrialOff(thisTrial), TZoneCol);0])==-1)]; % entry and exit indices w.r.t. Trial ON
        %thisTrialInZone = TrialInfo.Timestamps(thisTrial,1) + thisTrialInZone/SampleRate; % convert to seconds and offset w.r.t. trace start
        thisTrialInZone = thisTrialInZone/SampleRate; % convert to seconds
        if ~isempty(thisTrialInZone)
            TrialInfo.InZone(thisTrial) = { thisTrialInZone };
        else
            TrialInfo.InZone(thisTrial) = { [] };
        end
        
        %% Which Perturbation
        WhichPerturbation = mode( MyData(TrialOn(thisTrial):TrialOff(thisTrial), PerturbationCol(1)) );
        PerturbationValue = mode( MyData(TrialOn(thisTrial):TrialOff(thisTrial), PerturbationCol(2)) );
        if WhichPerturbation
            if WhichPerturbation < 5 % Fake target zone
                TrialInfo.Perturbation{thisTrial,1} = 'FakeZone';
                if isempty(find(TargetZones(:,2) == PerturbationValue))
                    TrialInfo.Perturbation{thisTrial,2} = PerturbationValue;
                else
                    TrialInfo.Perturbation{thisTrial,2} = find(TargetZones(:,2) == PerturbationValue);
                end  
            else
                switch WhichPerturbation
                    case 6   % feedback halt (old style) or feedback pause
                        if PerturbationValue == 5 % feedback halt (old style)
                            TrialInfo.Perturbation{thisTrial,1} = 'Halt-I';
                            if ~isempty(find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1))
                                HaltStart = 1;
                                HaltStop = find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1);
                                TrialInfo.Perturbation{thisTrial,2} = [HaltStart HaltStop];
                            end
                        end
                        if PerturbationValue == 6 % feedback pause
                            TrialInfo.Perturbation{thisTrial,1} = 'Halt-II';
                            if ~isempty(find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1))
                                HaltStart = find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==1);
                                HaltStop = find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1);
                                TrialInfo.Perturbation{thisTrial,2} = [HaltStart HaltStop];
                            end
                        end
                    case 300 % No Odor
                        TrialInfo.Perturbation{thisTrial,1} = 'NoOdor';
                    case 400 % flip map
                        TrialInfo.Perturbation{thisTrial,1} = 'FlipMap';
                    case 500 % location offset I
                        TrialInfo.Perturbation{thisTrial,1} = 'Offset-I';
                        TrialInfo.Perturbation{thisTrial,2} = PerturbationValue; % offset added
                    case {600, 700} % location offset II and III
                        if WhichPerturbation == 600
                            TrialInfo.Perturbation{thisTrial,1} = 'Offset-II';
                        else
                            TrialInfo.Perturbation{thisTrial,1} = 'Offset-III';
                        end
                        %TrialInfo.Perturbation(thisTrial,2) = {PerturbationValue}; % offset added
                        % get timestamps for offset start and feedback restart
                        % this is encoded in the InRewardZone Col - see GUI
                        if ~isempty(find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==1))
                            OffsetStart = ...
                                find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==1);
                            FeedbackStart = ...
                                find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1,1,'last');
                            % convert to seconds w.r.t. trial start
%                             OffsetStart = TrialInfo.PerturbationStart(thisTrial)/SampleRate;
%                             FeedbackStart = TrialInfo.FeedbackStart(thisTrial)/SampleRate;
                            OffsetStart = OffsetStart/SampleRate;
                            FeedbackStart = FeedbackStart/SampleRate;
                            
                            TrialInfo.Perturbation{thisTrial,2} = ...
                                [PerturbationValue OffsetStart FeedbackStart]; % offset added, offset start, feedback start w.r.t. trial start
                        end
                    case 800 % gain change
                        TrialInfo.Perturbation{thisTrial,1} = 'GainChange';
                        TrialInfo.Perturbation{thisTrial,2} = PerturbationValue; % new gain
                    case 1000 % halts (newer version - with location flip)
                        TrialInfo.Perturbation{thisTrial,1} = 'Halt-Flip';
                        if ~isempty(find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1))
                            HaltStart = find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==1);
                            HaltStop = find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1);
                            % HaltLocation = PerturbationValue
                            TrialInfo.Perturbation{thisTrial,2} = [HaltStart HaltStop PerturbationValue];
                        end
                    case 1510 % halts (newer version - with location flip) that also were used as open loop templates
                        TrialInfo.Perturbation{thisTrial,1} = 'Halt-Flip-Template';
                        if ~isempty(find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1))
                            HaltStart = find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==1);
                            HaltStop = find( diff([ MyData(TrialOn(thisTrial):TrialOff(thisTrial), RZoneCol); 0] )==-1);
                            % HaltLocation = PerturbationValue
                            TrialInfo.Perturbation{thisTrial,2} = [HaltStart HaltStop PerturbationValue];
                        end
                    case 1100 % block shift perturbations
                        TrialInfo.Perturbation{thisTrial,1} = 'BlockShift';
                        TrialInfo.Perturbation(thisTrial,2) = PerturbationValue; % shift amount
                    case 1400
                        if thisTrial>1
                            TrialInfo.Perturbation{thisTrial,1} = 'RuleReversal';
                            TrialInfo.Perturbation{thisTrial,2} = PerturbationValue - 1; %TFType during reversal
                        end
                    case 1500 % open loop template
                        TrialInfo.Perturbation{thisTrial,1} = 'OL-Template';
                    case 1600 % replay trial
                        TrialInfo.Perturbation{thisTrial,1} = 'OL-Replay';
                end
            end
        else
            TrialInfo.Perturbation{thisTrial,1} = [];
        end
    end
    
    LastTrialIdx = TrialOff(thisTrial); % current trial's end Idx
end

if any(strcmp(TrialInfo.Perturbation,'OL-Template'))
    %% hack : sometimes there are two templates but only one was replayed
    % lets assume its the second one for now
    % in such cases, there are no replay trials between the two templates
    TemplateTrials(:,1) = find(diff(strcmp(TrialInfo.Perturbation,'OL-Template'))== 1) + 1;
    if ~strcmp(TrialInfo.Perturbation(end,1),'OL-Template')
        TemplateTrials(:,2) = find(diff(strcmp(TrialInfo.Perturbation,'OL-Template'))==-1);
    else
        TemplateTrials(:,2) = [find(diff(strcmp(TrialInfo.Perturbation,'OL-Template'))==-1); size(TrialInfo.Perturbation,1)];
    end
    if size(TemplateTrials,1)>1
        first_replay = find(strcmp(TrialInfo.Perturbation,'OL-Replay'),1,'first');
        first_valid_template = find(TemplateTrials(:,2)<first_replay,1,'last');
        for t = 1:(first_valid_template-1)
            invalid_templates = TemplateTrials(t,1):TemplateTrials(t,2);
            for x = 1:numel(invalid_templates)
                TrialInfo.Perturbation{invalid_templates(x),1} = [];
            end
        end
    end
end
%% Extras: count trials of each target zone type
for i = 1:size(TargetZones,1)
    TargetZones(i,4) = numel( find(TrialInfo.TargetZoneType == i));
end

end