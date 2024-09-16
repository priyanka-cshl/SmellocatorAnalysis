function [] = MakePassiveSessionTraces(TuningMatFile)
    load(TuningMatFile);

    TrialVector = session_data.trace(:,find(strcmp(session_data.trace_legend,'trial_on')));
    TrialVector(TrialVector>0) = 1;
    
    ts = find(diff(TrialVector));
    if ~TrialVector(1)
        n = floor((length(ts)/2));
        idx = reshape(ts(1:2*n),2,n)';
        ts = session_data.timestamps(idx);

        % check that there was no clock drift
        if any(abs(TuningTTLs(:,2) - (ts(:,2) + TimestampAdjust.Passive))>0.04)
            disp('clock drift in ephys and behavior files');
            keyboard;
        end

    end

end