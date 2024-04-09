%% script to process a given tuning session where PID measurements were done
TuningFile = '/mnt/grid-hs/pgupta/Behavior/OBimaging/OBimaging_20220901_o5.mat'; % discrete tuning with PID measurements

%% basic processing of the tuning file
[MyTuningTrials, TuningSequence, ~, TuningSettings] = ParseTuningSession(TuningFile);
TuningSequence(1,:) = [];

% ugly hack to create the equivalent trial list using OEPS timestamps
TuningTTLs(:,1:3)   = MyTuningTrials(:,5:7); % TS-Trial-ON TS-Trial-OFF TS-Trial-duration
TuningTTLs(:,5)     = TuningSequence((0+ (1:size(TuningTTLs,1))),2); % odor ID
TuningTTLs(:,7)     = TuningSequence((0+ (1:size(TuningTTLs,1))),1); % odor location
TuningTTLs(:,8)     = (1:size(MyTuningTrials,1))'; % trial ID - fake it
TuningTTLs(:,12)    = (1:size(MyTuningTrials,1))'; % trial ID - fake it

[TuningInfo] = ParseTuningSubtrials(TuningTTLs, TuningSequence, TuningSettings);

%% extract the PID values/traces from the raw file
PIDChannel = 6; 

load(TuningFile);
PIDTrace(:,1) = session_data.timestamps;
PIDTrace(:,2) = session_data.trace(:,PIDChannel);

PIDvals = []; SteadyState = []; MeanVals = [];
for odor = 0:1:3
    
    whichperiods    = find(TuningInfo.Odor == odor);
    whichlocations  = TuningInfo.Location(whichperiods);
    whichtrials     = TuningInfo.TuningTrialID(whichperiods);
    for n = 1:numel(whichperiods)
        t = TuningInfo.trialtimes(whichperiods(n),1:2); % timestamps
        [~,idx1] = min(abs(PIDTrace(:,1)-t(1)));
        [~,idx2] = min(abs(PIDTrace(:,1)-t(2)));
        PIDvals{n,odor+1}  = PIDTrace(idx1:idx2,2);
        SteadyState(n,odor+1) = mean(PIDvals{n,odor+1}(end-99:end));
    end
    
    %%
    figure;
    count = 0;
    for loc = -90:15:90
        count = count + 1;
        subplot(1,13,count);
        hold on;
        f = find(whichlocations==loc);
        for i = 1:numel(f)
            plot(PIDvals{f(i),odor+1});
        end
        set(gca,'YLim',[-1.5 1]);
        MeanVals{odor+1}(count,1) = mean(SteadyState(f,odor+1));
        MeanVals{odor+1}(count,2) = std(SteadyState(f,odor+1));
    end
end
    