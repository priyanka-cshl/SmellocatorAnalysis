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
                    TrialVector(idx1:idx2) = -4; % passive tuning
%                     if ~ITIAirState
%                         % also flag time when manifold air actually turned on
%                         keyboard;
%                     end
                else % air trial
                    % no need to do anything to the OdorVector -
                    % it is already set to 4 as default;
                    [~,idx1] = min(abs(timestamps-TuningTTLs(t,1))); % trial start
                    [~,idx2] = min(abs(timestamps-TuningTTLs(t,2))); % trial end
                    TrialVector(idx1:idx2) = -4; % passive tuning
                end
            end
        else
            % get the TTL timings from replay TTLs
            whichReplay = find(ReplayTTLs.TrialID==TuningTTLs(t,8));
            ValveTS = ReplayTTLs.OdorValve{whichReplay}(:,1:2) + TuningTTLs(t,1);
            for n = 1:size(ValveTS,1) % every subtrial
                [~,idx1] = min(abs(timestamps-ValveTS(n,1))); % trial start
                [~,idx2] = min(abs(timestamps-ValveTS(n,2))); % trial end
                OdorVector(idx1:idx2) = ReplayTTLs.OdorValve{whichReplay}(n,4); % odor identity
                TrialVector(idx1:idx2) = -3; % replays
            end
        end
    end

end

PassiveOut.Timestamps{1} = timestamps;
PassiveOut.Odor{1} = OdorVector;
PassiveOut.Trial{1} = TrialVector;

%% tally valve timings with ephys TTLs
for x = 1:3
    odorvector = PassiveOut.Odor{1};
    odorvector(odorvector~=x) = 0;
    odorTS = PassiveOut.Timestamps{1}(find(diff(odorvector)));
    odorTS = reshape(odorTS,2,floor(numel(odorTS)/2))';

    [~,t1] = min(abs(TTLs.(['Odor',num2str(x)])(:,1)-odorTS(1,1)));  
    t2 = t1 + size(odorTS,1) - 1;

    if any(abs(odorTS(:,1)-TTLs.(['Odor',num2str(x)])(t1:t2,1))>0.005) || ...
       any(abs(odorTS(:,2)-TTLs.(['Odor',num2str(x)])(t1:t2,2))>0.005)  
        keyboard;
    end
end

% for the manifold air
manifoldVector = 0*PassiveOut.Odor{1};
manifoldVector(manifoldVector==4) = 0;
t1 = min(find(TTLs.AirManifold(:,1)>=timestamps(1),1,'first'), ...
    find(TTLs.AirManifold(:,2)>=timestamps(1),1,'first'));
t2 = size(TTLs.AirManifold,1);
ValveTS = TTLs.AirManifold(t1:t2,:);

for n = 1:size(ValveTS,1) % every transition
    [~,idx1] = min(abs(timestamps-ValveTS(n,1))); % start
    [~,idx2] = min(abs(timestamps-ValveTS(n,2))); % end
    manifoldVector(idx1:idx2) = 1;
end
PassiveOut.Manifold{1} = manifoldVector;

%% other traces
PassiveOut.Lever{1} = session_data.trace(:,find(strcmp(session_data.trace_legend,'lever_DAC')));
PassiveOut.Motor{1} = session_data.trace(:,find(strcmp(session_data.trace_legend,'stimulus_location_scaled')));
PassiveOut.Rewards{1} = session_data.trace(:,find(strcmp(session_data.trace_legend,'rewards')));
PassiveOut.Sniffs{1} = session_data.trace(:,find(strcmp(session_data.trace_legend,'thermistor')));

%% sniffs
% add a filtered sniff trace
PassiveOut.SniffsFiltered{1}     = FilterThermistor(PassiveOut.Sniffs{1});
% add a digital sniff trace
load(WhereSession,'CuratedPassiveSniffTimestamps');
if ~exist('CuratedPassiveSniffTimestamps','var')
    ProcessSniffTimeStamps_GUI(PassiveOut,WhereSession);
else
    if size(CuratedPassiveSniffTimestamps,2) < 10
        CuratedPassiveSniffTimestamps(:,10) = 0;
    end
    LocationSniffs = PassiveOut.SniffsFiltered{1}*nan;
    DigitalSniffs = PassiveOut.SniffsFiltered{1}*0;
    for n = 1:size(CuratedPassiveSniffTimestamps)
        idx = CuratedPassiveSniffTimestamps(n,8:9);
        if CuratedPassiveSniffTimestamps(n,10) == -1
            DigitalSniffs(idx(1):idx(2)) = -1;
        else
            DigitalSniffs(idx(1):idx(2)) = 1;
        end
        if ~any(isnan(CuratedPassiveSniffTimestamps(:,4)))
            location = CuratedPassiveSniffTimestamps(n,4);
            LocationSniffs(idx(1):idx(2)) = location;
        end
    end

    PassiveOut.SniffsDigitized{1} = DigitalSniffs;
    PassiveOut.SniffsLocationed{1} = LocationSniffs;
end

%% flagging transient passive perturbations

end