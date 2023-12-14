function [EphysTuningTrials_all] = AlignPassiveTuningTrials_Multi(TuningTrials_all, TTLs, BehaviorTrials, TrialSequence_all)
global FileOrder; % order in which behavior and Tuning sessions happened
ClosedLoops = find(FileOrder>0);
ClosedLoops(end+1) = numel(FileOrder)+1; % just a hack
EphysTuningTrials_all = [];

%% independently process each tuning session
U = unique(TuningTrials_all(:,end));
TuningSegments = [BehaviorTrials(:,2)+1 [BehaviorTrials(2:end,1)-1; size(TTLs.Trial,1)]] ;
for n = 1:size(TuningSegments,1)
    EphysTuningTrials = TTLs.Trial(TuningSegments(n,1):TuningSegments(n,2),:);
    EphysTuningTrials(:,8) = TuningSegments(n,1):TuningSegments(n,2);
    
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
    
    % match trials across tuning files from ephys and from behavior
    % which tuning files to worry about
    whichTuningSess = FileOrder(find(FileOrder(ClosedLoops(n):(ClosedLoops(n+1)-1))<0) + ClosedLoops(n) - 1);
    if numel(whichTuningSess)>1
        % more than one tuning session
        % lets assume that this is because something went wrong and there a
        % few extra trials that need to be accounted for
        % step 1: take the later session
        whichTuningSess(1:end-1,:) = [];
        TuningTrials = TuningTrials_all(find(TuningTrials_all(:,end)==abs(whichTuningSess)),1:end-1);
        TrialSequence = TrialSequence_all(find(TrialSequence_all(:,end)==abs(whichTuningSess)),1:end-1);
        % delete early trials in the Ephys side that may be extra
        [r,lags] = xcorr((EphysTuningTrials(:,3)),(TuningTrials(:,7)));
        todelete = lags(find(r==max(r)));
        % sanity check
        durationdiff = EphysTuningTrials(todelete+(1:2),3) - TuningTrials(1:2,7);
        if any(abs(durationdiff)>0.01)
            disp('unexplained session weirdness');
            keyboard;
        else
            EphysTuningTrials(1:todelete,:) = [];
        end
        EphysTuningTrials_all{abs(whichTuningSess)} = BasicPipeline(TrialSequence,TuningTrials,EphysTuningTrials);
    else
        TuningTrials = TuningTrials_all(find(TuningTrials_all(:,end)==abs(whichTuningSess)),1:end-1);
        TrialSequence = TrialSequence_all(find(TrialSequence_all(:,end)==abs(whichTuningSess)),1:end-1);
        EphysTuningTrials_all{abs(whichTuningSess)} = BasicPipeline(TrialSequence,TuningTrials,EphysTuningTrials);
    end
end

function [EphysTuningTrials] = BasicPipeline(TrialSequence,TuningTrials,EphysTuningTrials)
    % lets try the simplest thing first 
    [r,lags] = xcorr((EphysTuningTrials(:,3)),(TuningTrials(:,7)));
    extratrials = lags(find(r==max(r)));
    if numel(extratrials) == 1 && extratrials < (size(EphysTuningTrials,1) - size(TuningTrials,1))
        if ~any(abs(TuningTrials(1:end-1,7) - EphysTuningTrials(extratrials+(1:size(TuningTrials,1)-1),3))>0.01) % last trial seems to be shorter sometimes
            % delete the extra trial
            EphysTuningTrials(1:extratrials,:) = [];
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
                disp('odor sequences matched in ephys and behavior files');
                % delete any extra trials in the Ephys side
                EphysTuningTrials(size(TrialSequence,1)+1:end,:) = [];
                % copy over the motor locations from TrialSequence
                EphysTuningTrials(:,7) = TrialSequence(:,1);
            else
                disp('sequence mismatch in ephys and behavior tuning files');
                keyboard;
                % try deleting one trial on the ephys side
                if ~any([TrialSequence(:,2) - EphysTuningTrials(2:(size(TrialSequence,1)+1),5)])
                    EphysTuningTrials(1,:) = [];
                    disp('odor sequences now matched in ephys and behavior files');
                    % delete any extra trials in the Ephys side
                    EphysTuningTrials(size(TrialSequence,1)+1:end,:) = [];
                    % copy over the motor locations from TrialSequence
                    EphysTuningTrials(:,7) = TrialSequence(:,1);
                end
                
            end
        end
    else
        alignment = 0;
        while alignment<=0
            % Looks like if there were passive replays - spurious trials get inserted
            % 1. match trial durations across TuningTrials and EphysTuningTrials
            if ~any(abs(TuningTrials(1:end-1,7) - EphysTuningTrials(1:size(TuningTrials,1)-1,3))>0.01) % last trial seems to be shorter sometimes
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
                        disp('sequence mismatch in ephys and behavior tuning files');
                        keyboard;
                    end
                    % find the first problematic replay
                    f = find(isnan(EphysTuningTrials(:,5)));
                    xx = f(find(f<Mismatches(1),1,'last'));
                    % try deleting the trial after the replay
                    temp = EphysTuningTrials(:,5);
                    temp(xx+1,:) = [];
                    if numel(TuningListMismatches(TrialSequence(:,3),temp))<numel(Mismatches)
                        EphysTuningTrials(xx+1,:) = [];
                        TuningTrials(xx+1,:) = []; % added PG 2023/09/28
                        Mismatches = TuningListMismatches(TrialSequence(:,3),EphysTuningTrials(:,5));
                    end
                end
                disp('odor sequences matched in ephys and behavior tuning files');
                % delete any extra trials in the Ephys side
                EphysTuningTrials(size(TrialSequence,1)+1:end,:) = [];
                TuningTrials(size(TrialSequence,1)+1:end,:) = []; % added PG 2023/09/28
                % delete any extra trials in the Behavior side
                TrialSequence(size(EphysTuningTrials,1)+1:end,:) = [];
                % copy over the motor locations from TrialSequence
                EphysTuningTrials(:,7) = TrialSequence(:,1);
                alignment = 1;
            else
                if alignment == 0
                    % PG - 23/03/16: sometimes there are samples dropped in
                    % NI which causes misestimation of trial duration by
                    % exactly 50 ms
                    weirdo = find(abs(TuningTrials(:,7) - EphysTuningTrials(1:size(TuningTrials,1),3))>0.01);
                    if ~any(round((TuningTrials(weirdo,7) - EphysTuningTrials(weirdo,3)),2,'decimal')~=0.05)
                        TuningTrials(weirdo,7) = EphysTuningTrials(weirdo,3);
                    else
                        if numel(weirdo)>10 && weirdo(1) == 1 && size(EphysTuningTrials,1) > size(TuningTrials,1)
                            % sometimes there seem to be excess trials in
                            % the beginning
                            nn = size(EphysTuningTrials,1) - size(TuningTrials,1);
                            if ~any(abs(TuningTrials(:,7) - EphysTuningTrials(nn:end-1,3))>0.01)
                                disp('warning: using a recent hack that was written for S6_20230710');
                                keyboard;
                                EphysTuningTrials(1:nn-1,:) = [];
                                EphysTuningTrials(end,:) = [];
                            end
                        elseif any(find(TrialSequence(:,1)==800))
                            U1 = unique(round(TuningTrials(:,7),0,'decimal'));
                            U2 = unique(round(EphysTuningTrials(:,3),0,'decimal'));
                            culprit = U2(find(~ismember(U2,U1,'rows')));
                            spuriousTrials = find(round(EphysTuningTrials(:,3),0,'decimal') == culprit);
                            EphysTuningTrials(spuriousTrials,:) = [];
                        else
                            EphysTuningTrials(1,:) = [];
                        end
                        alignment = -1;
                    end
                else
                    disp('ephys file has extra trials!');
                    keyboard;
                end
            end
        end
    end
    EphysTuningTrials(:,5) = EphysTuningTrials(:,5) - 1; % make odor identities as 0-3 instead of 1-4
    % append three more columns, Trial start, Trial Stop, Trial duration -
    % from the behavior files - useful for sniff alignment etc. later
    EphysTuningTrials(:,9:11) = TuningTrials(:,5:7);
end

end
