function [] = SmellocatorHelp(QueryWord)

% split query word into subfields
foo = regexp(QueryWord, '\.', 'split');

switch foo{1}
    case 'TrialInfo'
        disp('Summary info struct created from the behavioral recordings in Matlab');
        switch foo{2}
            case 'TrialID'
                disp('[1 x N] : Original Trial ID from the raw behavior session');
            case 'Offset'
                disp('[1 x N] : offset in trial start computed from lever trace w.r.t. digitally recorded trial start.');
                disp('This happens due to unexpected sample drops only in the analog lines on certain matlab versions');
            case 'TraceIndices'
                disp('[N x 4] : start, stop idx of the trace w.r.t. raw session indices.');
                disp('Ignore col3,4 if there are no sample drops in analog lines w.r.t. digital lines');
                disp('If samples were dropped, col 1,2 are for digital traces; 3,4 are for analog traces.');
            case 'TraceDuration'
                disp('[N x 1] : trace length in seconds');
            case 'SessionIndices'
                disp('[N x 4] : start, stop idx of the trial w.r.t. raw session indices.');
                disp('Ignore col3,4 if there are no sample drops in analog lines w.r.t. digital lines');
                disp('If samples were dropped, col 1,2 are for digital traces; 3,4 are for analog traces.');
            case 'SessionTimestamps'
                disp('[N x 4] : start, stop timestamps of the trial w.r.t. raw session timestamps.');
                disp('Ignore col3,4 if there are no sample drops in analog lines w.r.t. digital lines');
                disp('If samples were dropped, col 1,2 are for digital traces; 3,4 are for analog traces.');
            case 'TimeIndices'
                disp('[N x 2] : Trial ON, OFF idx w.r.t. this trace start');
            case 'TimeStamps'
                disp('[N x 2] : Trial ON, OFF timestamps w.r.t. this trace start');
            case 'Duration'
                disp('[N x 1] : trial (not trace) length in seconds');
            case 'Odor'
                disp('[N x 1] : odor identity for this trial');
            case 'OdorStart'
                disp('[N x 2] : odor start in seconds w.r.t. trial start');
                disp('Col1 is computed from the InRewardZone channel - digital line');
                disp('Col2 is computed from the Lever trace - analog line');
                disp('Col1,2 should be the same if there are no Analog vs. digital sample drops');
            case 'TargetZoneType'
                disp('[N x 1] : which target zone out of the ones listed in TargetZones');
            case 'TransferFunctionLeft'
                disp('[N x 1] : 1 = trial started with +ve motor locations');
            case 'Reward'
                disp('{1 x N}[1, 2, ... n] : reward timestamps w.r.t. trace start until next trial');
                disp('rewards before trial start have -ve timestamps');
            case 'Success'
                disp('[N x 1] : (nRewards > 1) for this trial');
            case 'InZone'
                disp('{1 x N}[1, 2, ... n] : Target Zone entry and exit timestamps w.r.t. trial start, in seconds');
            case 'Perturbation'
                disp('{N x 1} : Perturbation type');
                disp('''Fake Zone'': which target zone was the fake zone');
                disp('''Halt-I'': [odor halt start, odor halt stop], odor halt start = trial start');
                disp('''Halt-II'': [odor halt start, odor halt stop], odor halt starts when lever crosses certain thresh');
                disp('''Halt-Flip'': [odor halt start, odor halt stop, halt location], odor halt starts when lever crosses certain thresh');
                disp('''NoOdor'': clear air instead of odor while Trial ON');
                disp('''FlipMap'': Odor reverses direction mid-trial, when lever crosses certain thresh');
                disp('''Offset-I'': [Offset added]');
                disp('              Odor is laterally offset after lever in TargetZone >0.5 * hold duration');
                disp('''Offset-II'': [Offset added, offset start, and feedback start w.r.t. trial start]');
                disp('               Odor is laterally offset and held in place after lever in TargetZone >0.5 * hold duration')
                disp('               Feedback is activated only when lever is moved past a certain threshold');
                disp('''Offset-III'': [Offset added, offset start, and feedback start w.r.t. trial start]');
                disp('                same as Offset-II, but with two different possible offsets')
                disp('''GainChange'': [New Gain], increase or decrease in lever-to-odor gain');
                disp('''BlockShift'': [shift size], rewarded location is shifted off-center');
                disp('''OL-Template'': close loop trial that was used as template for replay');
                disp('''OL-Replay'': replay trial');
            otherwise
                disp('fieldname unrecognized: check spelling or ask Priyanka');
        end
    case {'Traces','PassiveReplayTraces'}
        disp('Continuous traces for various behavioral variables for each trial, extracted from the matlab behavior/tuning files');
        disp('Traces span from ''startoffset'' seconds before trial start and go upto next trial start');
        disp('Two consecutive trial traces will thus overlap by ''startoffset'' seconds');
        switch foo{2}
            case 'Lever'
                disp('Recaled Lever trajectory (0 to 5 Volts), 5 = close to body');
            case 'Motor'
                disp('Odor Location. This is the calibrated position (-100 to 100) not the raw encoder signal.');
            case 'Sniffs'
                disp('Respiration signal - can be thermistor or pressure sensor');
            case 'Licks'
                disp('Lick detector state - read from the piezo comparator');
            case 'Trial'
                disp('Trial state. 0 = OFF, 1,2,3,4 = ON, value indicates which odor was ON, 4 is Air');
            case 'Rewards'
                disp('Water valve state');
            otherwise
                disp('fieldname unrecognized: check spelling or ask Priyanka');
        end
    case 'errorflags'
        disp('col1 : samples dropped on the analog channels, relative to the digital lines in matlab?');
        disp('col2 : timestamps dropped in general, in matlab?');
        disp('col3 : is Rotary encoder analog voltage to position calibration stable in this session?');
        disp(['col4 : was Motor position within ',char(177),'4 each time home sensor was HIGH?']);
        disp(['col5 : was Motor position within ',char(177),'10 each time InTargetZone indicator was HIGH?']);
    case 'TargetZones'
        disp('List of all target zones used in this session');
        disp('column 1-3: High Lim, Target, Low Lim in terms of Lever position');
        disp('column 4: number of trials per target zone type');
    case 'TuningTTLs'
        disp('Col1 : Trial ON Timestamp in OEPS timebase');
        disp('Col2 : Trial OFF Timestamp in OEPS timebase');
        disp('Col3 : Trial duration in seconds in OEPS timebase');
        disp('Col4 : Odor ON Timestamp in OEPS timebase');
        disp('Col5 : Odor Identity - 1 = Air, 2,3,4 are odors 1-3 respectively, NaN = passive replay trial');
        disp('Col6 : Odor OFF Timestamp in OEPS timebase');
        disp('Col7 : Motor location from behavior file, true location is within ',char(177),'5 units');
        disp('Col8 : Trial ID - w.r.t. all detetected Trials in OEPS and stored in TTLs');
    case 'TTLs'
        disp('Timestamps of digital events in OEPS file in OEPS timebase');
        disp('col1 : Timestamp of Event ON');
        disp('col2 : Timestamp of Event OFF');
        disp('col3 : Event duration, in seconds');
        switch foo{2}
            case 'Trial'
                disp('col4 : Time of Odor/Air ON w.r.t. Trial Start');
                disp('col5 : Odor Identity - 0 = Air, 1,2,3 are odors 1-3 respectively');
            otherwise
                disp('fieldname unrecognized: check spelling or ask Priyanka');
        end
    case 'ReplayTTLs'
        disp('Sequence of odor Valve transitions within each replay trial');
        disp('OdorValve{n}: Valve transitions for the nth replay w.r.t. replay trial start');
        disp('    Each row- ValveOn, ValveOff, Duration, Odor Identity'); 
        disp('TrialID{n} : Replay Trial ID in OEPS base - same as in the struct TTLs');
    case 'SingleUnits'
        disp('All putative single units from the corresponding recording session');
        switch foo{2}
            case 'id'
                disp('original cluster ID from phy');
            case 'tetrode'
                disp('which tetrode was this cluster on');
            case 'spikecount'
                disp('total spikes');
            case 'spikes'
                disp('spiketimes in seconds');
            case 'quality'
                disp('User defined cluster quality (from phy): 0 = noise, 1 = MUA, 2 = Good');
            case 'trialtags'
                disp('trial# of the trial (w.r.t. OEPS TTLs) that just started before this spike');
                disp('values are negative if it is a passive tuning trial (w.r.t. TuningTTLs');
            case 'trialalignedspikes'
                disp('spiketimes relative to the most recent trial''s'' start timestamp');
        end
    case 'startoffset'
        disp('Time window (in seconds) preceding trial start used for extracting ''Traces''');
    case 'SampleRate'
        disp('Behavior Sample Rate in Matlab (Hz)');
    case 'FileLocations'
        disp('Full File Paths for the relevant behavior and tuning file from Matlab, recording file from OEPS and sorted spikes from phy');
    otherwise
    	disp('Variable unrecognized: check spelling or ask Priyanka');              
end
        
end