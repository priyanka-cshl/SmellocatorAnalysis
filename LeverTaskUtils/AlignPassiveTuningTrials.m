function [EphysTuningTrials] = AlignPassiveTuningTrials(TuningTrials, TTLs, SkipTrials, TrialSequence)

if (size(TTLs.Trial,1) - SkipTrials) >= size(TuningTrials,1)
    EphysTuningTrials = TTLs.Trial(SkipTrials+1:end,1:3); % [TS-TrialOn TS-TrialOff Duration]
    % EphysTuningTrials = [TS-On TS-Off Duration ...
    %  TS-OdorON (w.r.t. TS-On) OdorID TS-OdorOFF ((w.r.t. TS-On) ...
    %  OdorLocation Trial# (from TTLs.Trial)]
    
    EphysTuningTrials(:,8) = (SkipTrials+1):size(TTLs.Trial,1);
    
    % first trial is crap - delete it
    EphysTuningTrials(1,:) = [];
    
    % Assign odor identities
    for i = 1:size(EphysTuningTrials,1)
        tstart = EphysTuningTrials(i,1);
        tstop = EphysTuningTrials(i,2);
        O1 = intersect(find(TTLs.Odor1(:,1)>tstart),find(TTLs.Odor1(:,1)<tstop));
        O2 = intersect(find(TTLs.Odor2(:,1)>tstart),find(TTLs.Odor2(:,1)<tstop));
        O3 = intersect(find(TTLs.Odor3(:,1)>tstart),find(TTLs.Odor3(:,1)<tstop));
        if ~numel(vertcat(O1,O2,O3))
            EphysTuningTrials(i,5) = 1;
        elseif numel(vertcat(O1,O2,O3)) == 1
            x = find([~isempty(O1) ~isempty(O2) ~isempty(O3)]); % which odor
            EphysTuningTrials(i,5) = 1 + x;
            EphysTuningTrials(i,[4 6]) = eval(['TTLs.Odor',num2str(x),'(O',num2str(x),',1:2);']);
        else
            EphysTuningTrials(i,5) = NaN; % passive replay
        end
    end
    
    % if there are no passive replays
    if ~any(find(TrialSequence(:,1)==999))
        % do motor locations in TuningTrials and TrialSequence match
        y = find(abs(TrialSequence(:,1) - TuningTrials(:,1))>5);
        if any(y)
            disp('location mismatches in .m tuning file')
            keyboard;
        else
            % do odor identities in EphysTuningTrials and TrialSequence match
            if ~any([TrialSequence(:,2) - EphysTuningTrials(1:size(TrialSequence,1),5)])
                disp('odor sequences match in ephys and behavior files');
                % delete any extra trials in the Ephys side
                EphysTuningTrials(size(TrialSequence,1)+1:end,:) = [];
                % copy over the motor locations from TrialSequence
                EphysTuningTrials(:,7) = TrialSequence(:,1);
            else
                disp('sequence mismatch in ephys and behavior tuning files');
                keyboard;
            end
        end
    else
        % Looks like if there were passive replays - spurious trials get inserted
        % 1. match trial durations across TuningTrials and EphysTuningTrials
        if ~any(abs(TuningTrials(:,7) - EphysTuningTrials(1:size(TuningTrials,1),3))>0.005)
            % delete any extra trials in the Ephys side
            EphysTuningTrials(size(TuningTrials,1)+1:end,:) = [];
            % copy over the motor locations from TuningTrials
            EphysTuningTrials(:,7) = TuningTrials(:,1);
            % delete spurious trials - if any
            % to ignore any replay trials and the one just after - convert
            % odor identities in the behavior file to NaN
            TrialSequence(:,3) = TrialSequence(:,2);
            TrialSequence(TrialSequence(:,1)==999,3) = NaN;
            TrialSequence(1+find(TrialSequence(:,1)==999),3) = NaN;
            % look for odor identity mismatches that are not NaNs
            Mismatches = TuningListMismatches(TrialSequence(:,3),EphysTuningTrials(:,5));
            while ~isempty(Mismatches) 
                % make sure there are extra trials in the ephys list that
                % need to be deleted
                if size(EphysTuningTrials,1)<size(TrialSequence,1)
                    disp('sequence mismatch in ephys and behavior tunig files');
                    keyboard;
                end
                % find the first problematic replay
                f = find(isnan(EphysTuningTrials(:,5)));
                x = f(find(f<Mismatches(1),1,'last'));
                % try deleting the trial after the replay 
                temp = EphysTuningTrials(:,5);
                temp(x+1,:) = [];
                if numel(TuningListMismatches(TrialSequence(:,3),temp))<numel(Mismatches)
                    EphysTuningTrials(x+1,:) = [];
                    Mismatches = TuningListMismatches(TrialSequence(:,3),EphysTuningTrials(:,5));
                end
            end
            disp('odor sequences match in ephys and behavior tuning files');
            % delete any extra trials in the Ephys side
            EphysTuningTrials(size(TrialSequence,1)+1:end,:) = [];
            % copy over the motor locations from TrialSequence
            EphysTuningTrials(:,7) = TrialSequence(:,1);
        else
            disp('ephys file has extra trials!');
            keyboard;
        end
    end
    EphysTuningTrials(:,5) = EphysTuningTrials(:,5) - 1; % make odor identities as 0-3 instead of 1-4
else
    EphysTuningTrials = [];
end
