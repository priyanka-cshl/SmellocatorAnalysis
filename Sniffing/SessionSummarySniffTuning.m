MySession = '/mnt/data/Processed/Behavior/Q4/Q4_20221109_r0_processed.mat';

%% get the data loaded
[TracesOut, ColNames, TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

%% Get all spikes, all units aligned to trials
[AlignedSniffs, sniffAlignedSpikes, trialAlignedSpikes, whichtetrode] = ...
    SniffAlignedSpikeTimes(SingleUnits,TTLs,size(TrialInfo.TrialID,2),TrialInfo,MySession);

%% Loop through units & odors
for whichunit = 1:size(SingleUnits,2)
    for whichodor = 1:3
        [~,FR{whichunit,whichodor},BinOffset] = PlotSortedSniffs(whichunit, whichodor, trialAlignedSpikes, AlignedSniffs, ...
            TrialInfo, 'plotevents', 0, 'plotspikes', 0, 'psth', 1, ...
            'sortorder', 0, 'alignto', 1, 'warptype', 0);
    end
end

%% make a dumb matrix
for whichunit = 1:size(SingleUnits,2)
    whichodor = 1;
    snifftype = 1; % ITI
    foo = FR{whichunit,whichodor}{1};
    % delete the previous sniff points
    todelete = (abs(BinOffset)/2) - 100;
    foo(1:todelete,:) = [];
    
        % zscore 
    foo = (foo - mean(foo))/std(foo);
    
    [a,b] = max((foo(101:300)));
    Metrics(whichunit,2) = b + 100;
    Metrics(whichunit,1) = foo(b);
    AirTuning(whichunit,:) = foo/a;

end

% sort by response amplitude and latency
[~,unitorder] = sortrows(Metrics,[2 1]);


%plot((1:size(FR{t},1))*0.002+BinOffset/1000,FR{t},'Linewidth',2,'Color','k');