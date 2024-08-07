function [MyData] = PatchSampleDrops(MyData, myEphysdir)
%% function to patch samples in the analog traces that got dropped
% from the parallel OEPS recordings

%% read the open ephys flat binary file
session = Session(myEphysdir);
[TTLs] = OEPS_getTTLs_newGUI(session);

TS = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps;
TS = TS - session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps(1);

% first find the timestamp adjustment that needs to be made
TrialOn.behavior = MyData(find(diff(MyData(:,6))>0,5,'first'),1);

if ~any((diff(TrialOn.behavior) - diff(TTLs.Trial(1:5,1))) > 0.05)
    TimestampAdjust.ClosedLoop = median(TTLs.Trial(1:5,1) - TrialOn.behavior);
    % make sure that OEPS didn't crash midway
    if ~any(TS>MyData(end,1)+TimestampAdjust.ClosedLoop)
        disp('OEPS crashed mid-session - aborting patch sample drops');
        return;
    end
else
    disp('could not match trial times to calculate timestampadjust');
    keyboard;
end

%%
drop_point = find(abs(diff(MyData(:,1)))>0.01, 1, 'first'); % first indices at which timestamps were dropped
while ~isempty(drop_point)
    [MyData] = patch_smellocatordata_gap(drop_point,MyData,session,TS,TimestampAdjust.ClosedLoop,TTLs);
    drop_point = find(abs(diff(MyData(:,1)))>0.01, 1, 'first'); % first indices at which timestamps were dropped
end

%LeverSnippet


    function [MyData] = patch_smellocatordata_gap(drop_point,MyData,session,TS,TSadjuster,TTLs)
        % drop_points = last index in MyData after which the samples got dropped.
        % MyData = actual old data with the timestamp drops
        
        % 1. create a dummy block
        ts_before   = MyData(drop_point,1);
        ts_after    = MyData(drop_point+1,1);
        ts_dropped  = (ts_before + 0.002):0.002:(ts_after - 0.002);
        DummyBlock  = nan(numel(ts_dropped),size(MyData,2));
        
        % 2. find trial off indices and timestamps
        idx_trialoff    = find(MyData((drop_point+1):end,6)==0,1,'first') + drop_point;
        ts_full         = [ts_dropped MyData((drop_point+1):idx_trialoff,1)'];
        
        % 3. convert to OEPS indices to get the correct trace portions out
        % recalculate TSadjuster
        [timediff,whichtrial] = min(abs(TTLs.Trial(:,2)-(ts_full(end)+TSadjuster)));
        if timediff > 0.0025
            disp('correcting for clock offset drift..')
            TSadjuster = TTLs.Trial(whichtrial,2)-ts_full(end);
        end
        
        [~,idx_OEPS(1)] = min(abs(TS-ts_full(1)-TSadjuster));
        [~,idx_OEPS(2)] = min(abs(TS-ts_full(end)-TSadjuster));
        idx_ephys = idx_OEPS(1):idx_OEPS(2);
        ts_OEPS   = TS(idx_ephys);
        % Lever, RE, Thermistor, Piezo
        OEPSchanIDs = [43 45 41 44];
        Traces_OEPS_raw = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(OEPSchanIDs,idx_ephys)';
        Traces_OEPS = []; 
        for ch = 1:numel(OEPSchanIDs)
            % convert to volts
            VoltMultiplier = session.recordNodes{1}.recordings{1}.info.continuous.channels(OEPSchanIDs(ch)).bit_volts;
            Traces_OEPS(:,ch) = double(Traces_OEPS_raw(:,ch))*VoltMultiplier;
        end
        %Traces_OEPS = 2.5 + (double(Traces_OEPS)/12800);
        
        % 4. interpolate to get traces at behavior resolution
        % Vq = interp1(X,V,Xq)
        Traces_interp = interp1( (ts_OEPS-ts_OEPS(1)), Traces_OEPS, (ts_full-ts_full(1)) );
        
        % 5. get behavior traces to confirm the match
        whichcolumns = [4 5 15 18]; % Lever, RE, Thermistor, Piezo
        trial_column = 6;
        In_RZ_column = 8;
        Data_B = MyData((drop_point+1):idx_trialoff,whichcolumns);
        % which portion of the interpolated traces to compare?
        [~,idx_split] = min(abs(ts_full-ts_after));
        Data_E = Traces_interp(idx_split:end,:);
        
        % 6. calculate the DC offset
        trace_offset = nanmedian(Data_E-Data_B,1);
        
        % 7. compare
        for n = 1:4
            if nanmedian(abs(Data_E(:,n)-trace_offset(n) - Data_B(:,n))) < 0.025
                if n == 3
                    disp('thermistor matched!');
                end
                % good match
                [~,~,which_idx] = intersect(ts_dropped,ts_full);
                DummyBlock(:,whichcolumns(n)) = Traces_interp(which_idx,n) - trace_offset(n);
                if any(isnan(DummyBlock(:,1)))
                    % add timestamps
                    DummyBlock(:,1) = ts_dropped;
                    
                    % fix trial ON vector
                    trial_idx_OEPS      = find(TTLs.Trial(:,1)>=ts_OEPS(1),1,'first');
                    TrialOn_OEPS        = TTLs.Trial(trial_idx_OEPS,1);
                    TrialOn_Behavior    = TrialOn_OEPS - TSadjuster;
                    trial_idx_Behavior  = find(ts_full>=TrialOn_Behavior,1,'first');
                    DummyBlock(:,trial_column) = MyData(drop_point,trial_column);
                    DummyBlock(trial_idx_Behavior:end,trial_column) = MyData(drop_point+1,trial_column);
                    
                    % fix InRewardZone vector
                    Manifold_idx_OEPS   = find(TTLs.AirManifold(:,1)<=TrialOn_OEPS,1,'last');
                    ManifoldOn_OEPS     = TTLs.AirManifold(Manifold_idx_OEPS,1);
                    ManifoldOn_Behavior = ManifoldOn_OEPS - TSadjuster;
                    Manifold_idx_Behavior = find(MyData(1:(drop_point + size(DummyBlock,1)),1)>=ManifoldOn_Behavior,1,'first') - drop_point;
                    if Manifold_idx_Behavior>0
                        if MyData(drop_point,In_RZ_column)
                            DummyBlock(:,In_RZ_column) = 1;
                            DummyBlock(Manifold_idx_Behavior:end,In_RZ_column) = 0;
                        elseif Manifold_idx_Behavior == 1 && MyData(drop_point-1,In_RZ_column)
                            DummyBlock(:,In_RZ_column) = 1;
                            DummyBlock(Manifold_idx_Behavior:end,In_RZ_column) = 0;
                        else
                            keyboard;
                        end
                    elseif ~MyData(drop_point,In_RZ_column)
                        DummyBlock(:,In_RZ_column) = 0;
                    else
                        keyboard;
                    end
                    
                    % fix InTargetZone(7), Rewards and(9) HomeSensor cols (14)
                    % just put them all to zero
                    % Lever can't be in target zone post init and pretrial start
                    % no rewards in this window
                    % home not possible either
                    DummyBlock(:,[7 9 14]) = 0;
                    
                    % Missing Licks - patch in from OEPS
                    missinglicksOn = intersect(find(TTLs.Licks(:,1)>=ts_dropped(1)+TSadjuster), ...
                        find(TTLs.Licks(:,1)<=ts_dropped(end)+TSadjuster));
                    missinglicksOff = intersect(find(TTLs.Licks(:,2)>=ts_dropped(1)+TSadjuster), ...
                        find(TTLs.Licks(:,2)<=ts_dropped(end)+TSadjuster));
                    if isempty(missinglicksOn) && isempty(missinglicksOff)
                        % no licks found
                        DummyBlock(:,10) = 0;
                    else
                        % licks found
                        AllLicks = vertcat([missinglicksOn ones(numel(missinglicksOn),1)],[missinglicksOff zeros(numel(missinglicksOff),1)]);
                        AllLicks = sortrows(AllLicks,1);
                        % are these legit detections
                        if ~AllLicks(1,2)==MyData(drop_point,10) %&& AllLicks(end,2)==MyData(drop_point+1,10)
                            for l = 1:size(AllLicks,1)
                                % get indices in the dropped trace
                                if AllLicks(l,2)
                                    [~,AllLicks(l,3)] = min(abs(ts_dropped+TSadjuster - TTLs.Licks(AllLicks(l),1)));
                                else
                                    [~,AllLicks(l,3)] = min(abs(ts_dropped+TSadjuster - TTLs.Licks(AllLicks(l),2)));
                                end
                            end
                            if mod(size(AllLicks,1),2)
                                AllLicks(end+1,:) = [NaN ~AllLicks(end,2) size(DummyBlock,1)];
                            end
                            for l = 2:2:size(AllLicks,1)
                                DummyBlock(:,10) = ~AllLicks(1,2);
                                DummyBlock(AllLicks(l-1,3):AllLicks(l,3),10) = AllLicks(1,2);
                            end
                        else
                            disp('lick transitions do not make sense');
                            % keyboard;
                            DummyBlock(:,10) = 0;
                        end
                                
                    end
                end
            else
                if (n == 1) || (n == 2)
                    disp('Lever/RE traces dont match up across the two files');
                    keyboard;
                elseif  n == 3 && ~any(isnan(DummyBlock(:,1)))
                    disp('thermistor mismatch - interpolating ...');
                    % thermistor trace didn't match up - lets interpolate
                    sniff_interp = linspace(MyData(drop_point,whichcolumns(n)),...
                        MyData(drop_point+1,whichcolumns(n)),...
                        size(DummyBlock,1)+2);
                    DummyBlock(:,whichcolumns(n)) = sniff_interp(2:end-1);
                else
                    disp('piezo mismatch - ignoring ...');
                end
            end
        end
        
        % 8. pad the block
        MyData = vertcat(MyData(1:drop_point,:), DummyBlock, MyData((drop_point+1):end,:));
        
    end

end
