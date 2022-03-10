%% for O3
SessionPath = 'O3/O3_20210927_r0_processed.mat';

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

%% Concatenate traces and get one matrix with all behavior variables
% SampleRate = behavior sample rate;
[TracesOut, ColNames] = ConcatenateTraces2Mat(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);

% pad the samples missing from raw behavior file to prevent indexing mismatches
TracesOut = vertcat(zeros(TrialInfo.TraceIndices(1,1),size(TracesOut,2)), TracesOut);

%% sniff trace
%sniffs = TracesOut(:,3);
fband = [0.5 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
TracesOut(:,3) = filtfilt(b,a,TracesOut(:,3)); %apply the filter to x(t)

% find peaks
[pks,pkloc] = findpeaks(TracesOut(:,3),'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);
% find valleys
[vls,vlloc] = findpeaks(-TracesOut(:,3),'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);

% sanity checks
for i = 1:numel(pks)-1 
    % is there a unique valley
    thispkvalley = intersect(find(vlloc>pkloc(i)),find(vlloc<pkloc(i+1)));
    if numel(thispkvalley==1)
        Loc(i,:) = [pkloc(i) vlloc(thispkvalley)];
    else
        keyboard;
    end
end

%%
figure; 
subplot(1,3,1); 
histogram(diff(Loc(:,1))); 
subplot(1,3,2); 
histogram(diff(Loc(:,2))); 
subplot(1,3,3); 
histogram((Loc(:,2)-Loc(:,1))); 

%% Trial Aligned Sniff times?
[AlignedSniffs] =  TrialAlignedSniffTimes(Loc,TracesOut(:,5),SampleRate);

%% Trial Aligned spikeTimes
TrialInfo.TargetEntry = NaN*TrialInfo.Odor;
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

%% Plot
figure;
AlignTo = 2; 
switch AlignTo
    case {1,2,6}
        myXlim = [-1.2 6];
    case {3,4}
        myXlim = [-5.2 1];
    case 5
        myXlim = [-2.2 5];
end

whichUnit = 15; whichOdor = 1;
subplot(3,2,[1 3]); hold on
[FRs, BinOffset, whichTZ] = ... 
PlotHaltFlips(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo);
set(gca,'XLim',myXlim);

subplot(3,2,[5]); hold on
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);

whichUnit = 1;
subplot(3,2,[2 4]); hold on
[SRs, BinOffset, whichTZ] = ... 
PlotHaltFlips(whichUnit, whichOdor, AlignedSniffs, Events, TrialInfo, AlignTo);
set(gca,'XLim',myXlim);

subplot(3,2,[6]); hold on
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);
