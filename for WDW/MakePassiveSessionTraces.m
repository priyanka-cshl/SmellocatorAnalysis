function [PassiveOut] = MakePassiveSessionTraces(WhereSession)

%% Load the relevant variables for the clsed loop processed file
load(WhereSession, 'FileLocations', 'TuningTTLs', 'TTLs', 'ReplayTTLs');
% load the actual tuning .mat file
load(FileLocations.Tuning); % loads session_data

%     'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
%     'startoffset', 'errorflags', 'SampleRate', ...
%     'TTLs', 'ReplayTTLs',  'SingleUnits', ...
%     'TimestampAdjust');

%% Get trial On-Off times from the trial state vector
TrialVector = session_data.trace(:,find(strcmp(session_data.trace_legend,'trial_on')));
TrialVector(TrialVector>0) = 1;

ts = find(diff(TrialVector));
if ~TrialVector(1)
    n = floor((length(ts)/2));
    idx = reshape(ts(1:2*n),2,n)';
    ts = session_data.timestamps(idx);
    ts(:,3) = ts(:,2)-ts(:,1);
else
    keyboard;
end

% quick sanity check
if size(TuningTTLs,1) ~= size(ts,1)
    keyboard;
end

% calculate the required timestamp adjustment
ephys       = TuningTTLs(:,2);
behavior    = ts(:,2);
myfit = fit(behavior,ephys-behavior,'poly1');

TimestampAdjust.Passive(2) = myfit.p2;
TimestampAdjust.Passive(1) = myfit.p1;

% make an adjusted timestamps vector
timestamps = session_data.timestamps + session_data.timestamps*TimestampAdjust.Passive(1) + TimestampAdjust.Passive(2);

%% construct the Trial and Odor vectors
% was the air On or Off during ITI?
f = find(strcmp(session_data.legends,'ITIAirState'));
if ~isempty(f)
    ITIAirState = session_data.params(f); % 1 if air was on throughout ITI
    % manifold air only turns off in the ITIs between the replay trials
    % find the time when the manifold turned on before session start
    f = find(TTLs.AirManifold(:,2)<TuningTTLs(1,1),1,'last');
    if (TTLs.AirManifold(f+1,1)-TuningTTLs(1,1))<0 && (TTLs.AirManifold(f+1,2)-TuningTTLs(1,1))>0
        ManifoldOnTS = TTLs.AirManifold(f+1,1);
    else
        keyboard;
    end
    % time of passive session start =
    % find(TTLs.AirManifold(:,1)<TuningTTLs(1,1),1,'last');
else
    keyboard;
    ITIAirState = 0;
end

OdorVector = (timestamps*0) + (ITIAirState*4);
TrialVector = timestamps*0;

if ITIAirState
    % fixing trial type for trial 1
    if isnan(TuningTTLs(1,7))
        if any(TuningTTLs(:,7)==800)
            if abs(TuningTTLs(1,3) - median(TuningTTLs(find(TuningTTLs(:,7)==800),3))) < 0.01
                TuningTTLs(1,7) = 800;
            else
                keyboard;
            end
        else
            keyboard;
        end
    end

    for t = 1:size(TuningTTLs,1) % every trial
        if TuningTTLs(t,7) < 990 % not a replay
            if TuningTTLs(t,7) == 800
                if TuningTTLs(t,5) % odor ON
                    [~,idx1] = min(abs(timestamps-TuningTTLs(t,4))); % odor start
                    [~,idx2] = min(abs(timestamps-TuningTTLs(t,6))); % odor end
                    OdorVector(idx1:idx2) = TuningTTLs(t,5); % odor identity
                    TrialVector(idx1:idx2) = -4;
%                     if ~ITIAirState
%                         % also flag time when manifold air actually turned on
%                         keyboard;
%                     end
                else % air trial
                    % no need to do anything to the OdorVector -
                    % it is already set to 4 as default;
                    [~,idx1] = min(abs(timestamps-TuningTTLs(t,1))); % trial start
                    [~,idx2] = min(abs(timestamps-TuningTTLs(t,2))); % trial end
                    TrialVector(idx1:idx2) = -4;
                end
            end
        else
            % get the TTL timings from replay TTLs
            keyboard;
            whichReplay = find(ReplayTTLs.TrialID==TuningTTLs(t,8));
            % manifold is turning off 1 second after trial off! :(
    
        end
    end

end


% make a timestamp trace
PassiveOut.Timestamps{1} = timestamps;
%PassiveOut.Lever =
%     whichTraces{1} = 'Lever'; whichTraces{2} = 'Motor'; whichTraces{3} = 'Sniffs';
%     whichTraces{4} = 'Trial'; whichTraces{5} = 'Odor';
%     whichTraces{6} = 'Rewards'; whichTraces{7} = 'Timestamps';


% check that there was no clock drift
if any(abs(TuningTTLs(:,2) - (ts(:,2) + TimestampAdjust.Passive))>0.04)
    disp('clock drift in ephys and behavior files');
    keyboard;
end


end