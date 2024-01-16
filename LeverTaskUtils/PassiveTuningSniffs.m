function [SniffTS_passive] = PassiveTuningSniffs(TuningInfo,MySession)

% determine passive tuning type
if size(TuningInfo.extras.sequence,1)>2 && size(TuningInfo.extras.sequence,2) == 2 % old style passive tuning
    tuningtype = 1;
elseif any(TuningInfo.extras.sessionsettings(:,1)==800)
    tuningtype = 2; % pseudorandom tuning
end

% determine if Air was ON or OFF during ITI
if length(TuningInfo.extras.sessionsettings) == 12
    ITIAirState = TuningInfo.extras.sessionsettings(12);
else
    ITIAirState = 0;
end

% load the passive sniff timestamps and convert to OEPS timebase
load(MySession,'SniffTS_passive'); %, 'TimestampAdjust'); % sniff timestamps in behavior timebase
loaded_adjuster = 0;
while ~loaded_adjuster
    lastwarn('', ''); % reset the lastwarn message and id
    load(MySession,'TimestampAdjust'); % this might throw a warning because TimestampAdjust may not exist
    [~, warnId] = lastwarn(); % if a warning was raised, warnId will not be empty.
    if(isempty(warnId))
        Passive_Timestamp_adjust = TimestampAdjust.Passive;
        loaded_adjuster = 1;
    else
        % reprocess the session in question
        disp('Run preprocess again');
        keyboard;
        % continue if you want to quickly reprocess
        [~,F,ext] = fileparts(MySession);
        PreprocessSmellocatorData(strrep(F,'_processed',ext),1);
    end
end

SniffTS_passive(:,1:3) = SniffTS_passive(:,1:3) + Passive_Timestamp_adjust; % Sniff Times in OEPS timebase

% keep track of previous and next sniffs (for plotting)
SniffTS_passive(:,6) = SniffTS_passive(:,4); % locations
SniffTS_passive(2:end,7:10) = SniffTS_passive(1:end-1,1:4); % previous sniff
SniffTS_passive(1:end-1,11:14) = SniffTS_passive(2:end,1:4); % next sniff
SniffTS_passive(:,4:5) = NaN*SniffTS_passive(:,4:5);
SniffTS_passive(:,15) = SniffTS_passive(:,3) - SniffTS_passive(:,1); % sniff duration
SniffTS_passive(:,16) = SniffTS_passive(:,2) - SniffTS_passive(:,1); % inh duration
%SniffTS_passive(:,17) = 1:size(SniffTS_passive,1); % sniff ID

if ~isempty(SniffTS_passive)
    for whichtrial = 1:size(TuningInfo.extras.sequence,1)
        if (TuningInfo.extras.sequence(whichtrial,1) == 998) || (TuningInfo.extras.sequence(whichtrial,1) == 999) % replay trials
            replaystart = TuningInfo.TTLs(whichtrial,1);
            replaystop  = TuningInfo.TTLs(whichtrial,2);
            first_inhalation = find(SniffTS_passive(:,1)>=replaystart,1,'first'); % after replay start
            last_inhalation  = find(SniffTS_passive(:,2)<=replaystop ,1,'last' ); % before replay end
            SniffTS_passive(first_inhalation:last_inhalation,5) = -2;
            SniffTS_passive(first_inhalation:last_inhalation,17) = whichtrial;
            
        else % tuning trials
            if ~isnan(TuningInfo.extras.sequence(whichtrial,1))
                if whichtrial ~= 1
                    prevtrialoff = TuningInfo.TTLs(whichtrial-1,2);
                else
                    prevtrialoff = 0;
                end
                trialstart  = TuningInfo.TTLs(whichtrial,1);
                trialstop   = TuningInfo.TTLs(whichtrial,2);
                odorstart   = TuningInfo.TTLs(whichtrial,4);
                odorstop    = TuningInfo.TTLs(whichtrial,6);
                whichodor   = TuningInfo.TTLs(whichtrial,5);
                
                first_inhalation = find(SniffTS_passive(:,1)>=prevtrialoff,1,'first'); % after prev trial end
                last_inhalation  = find(SniffTS_passive(:,2)<=trialstop ,1,'last' ); % before this trial end
                SniffTS_passive(first_inhalation:last_inhalation,17) = whichtrial;
                
                last_inhalation  = find(SniffTS_passive(:,2)<=trialstart ,1,'last' ); % before this trial start
                if ~ITIAirState
                    SniffTS_passive(first_inhalation:last_inhalation,5) = -1; % ITI sniffs
                else
                    SniffTS_passive(first_inhalation:last_inhalation,5) = 0; % ITI sniffs == Air Sniffs
                end
                
                first_inhalation = find(SniffTS_passive(:,1)>=trialstart,1,'first'); % after this trial start
                last_inhalation  = find(SniffTS_passive(:,2)<=trialstop ,1,'last' ); % before this trial end
                SniffTS_passive(first_inhalation:last_inhalation,5) = 0; % air sniffs
                
                first_inhalation = find(SniffTS_passive(:,1)>=odorstart,1,'first'); % after odor start
                last_inhalation  = find(SniffTS_passive(:,2)<=odorstop ,1,'last' ); % before odor end
                SniffTS_passive(first_inhalation:last_inhalation,5) = whichodor;
            end
        end
    end
end

end