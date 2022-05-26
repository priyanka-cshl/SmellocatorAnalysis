%LeverTaskPhotometryAnalysisTest
clc;clear all;

whatmouse = 'DH2';

switch whatmouse

    case 'DH2'
        WhichSession = '2020-03-12_17-37-26';
        SessionPath = 'C:\Users\Marie\Documents\data\Data_Photometry_Diego\Lever_task\DH2';
        myKsDir = fullfile(SessionPath,WhichSession);

    case 'DH4' % needs debugging 
        WhichSession = '2020-03-10_14-18-48';
        SessionPath = 'C:\Users\Marie\Documents\data\Data_Photometry_Diego\Lever_task\DH4';
        myKsDir = fullfile(SessionPath,WhichSession);
end

%% Get Trial Timestamps from the OpenEphys Events file
filename = fullfile(myKsDir,'all_channels.events');

%% heck to avoid going through empty files
temp = dir(filename);
if ~temp.bytes
    TTLs = [];
    EphysTuningTrials = [];
    AuxData = [];
    disp('empty events file');
    return
end

[data, timestamps, info] = load_open_ephys_data(filename); % data has channel IDs
% for some reason there is something wrong about that timestamps

% adjust for clock offset between open ephys and kilosort
[offset] = AdjustClockOffset(myKsDir);
% offset = offset/OepsSampleRate;
timestamps = timestamps - offset;


% Get various events
TTLTypes = unique(data);
Tags = {'Air', 'Odor1', 'Odor2', 'Odor3', 'Trial', 'Reward', 'AirManifold', 'Licks'};
for i = 1:numel(TTLTypes)
    On = timestamps(intersect(find(info.eventId),find(data==TTLTypes(i))));
    Off = timestamps(intersect(find(~info.eventId),find(data==TTLTypes(i))));
    % delete the first off value, if it preceeds the On
    Off(Off<On(1)) = [];
    On(On>Off(end)) = [];
    
    if length(On)>length(Off)
        keyboard;
        foo = [On(1:end-1) Off]';
        goo = foo(:);
        On(find(diff(goo)<0,1,'first')/2,:) = [];
    end
        
    temp = [On Off Off-On];
    
    % ignore any transitions faster than 1 ms - behavior resolution is 2 ms
    temp(temp(:,3)<0.001,:) = [];
    TTLs.(char(Tags(i))) = temp;
end

%% Read Photoreceivers and LEDs data from open ephys file
% use what is written below
%     foo = dir(fullfile(myKsDir,'*_ADC1.continuous')); % pressure sensor
%     filename = fullfile(myKsDir,foo.name);
%     [Auxdata1, timestamps, ~] = load_open_ephys_data(filename); % data has channel IDs
%     foo = dir(fullfile(myKsDir,'*_ADC2.continuous')); % thermistor
%     filename = fullfile(myKsDir,foo.name);
%     [Auxdata2, ~, ~] = load_open_ephys_data(filename); % data has channel IDs

[modData_1, Timestamps_PR_1, info_PR_1]  = load_open_ephys_data('100_ADC3.continuous');
[modLED_1, Timestamps_LED_1, info_LED_1] = load_open_ephys_data('100_ADC5.continuous');

%% Constants

modFreq_1          = 211;       
modAmp_1           = 0.6;
samplingRate       = info_PR_1.header.sampleRate;
lowCutoff          = 15;

%% Prepare reference data and generate 90deg shifted reference data
   
shift_modLED_1     = modLED_1 - mean(modLED_1);                            % Remove DC offset
samplesPerPeriod   = (samplingRate/modFreq_1);
quarterPeriod      = round(samplesPerPeriod/4);
shift_modLED90_1   = circshift(shift_modLED_1,[1 quarterPeriod]);

%% Quadrature decoding and filtering                                       % Element-by-element array multiplication 
   
processedData0_1    = modData_1 .* shift_modLED_1;                         % 0 degrees data correction                          
processedData90_1  = modData_1 .* shift_modLED90_1;                        % 90 degrees data correction 

%% Low pass filter
    
norm_lowCutoff     = lowCutoff/(samplingRate/2);                           % CutOff normalized by half sampling rate 
[b, a]             = butter(5, norm_lowCutoff,'low');                      % '5th order' butterworth low pass filter

paddedData0_1        = processedData0_1(1:samplingRate,1);             
paddedData90_1       = processedData90_1(1:samplingRate,1);
demodDataFilt0_1     = filtfilt(b,a,[paddedData0_1; processedData0_1]);    % pad the data to suppress windows effect upon filtering
demodDataFilt90_1    = filtfilt(b,a,[paddedData90_1; processedData90_1]);        
processedData0f_1    = demodDataFilt0_1(samplingRate + 1: end, 1);
processedData90f_1   = demodDataFilt90_1(samplingRate + 1: end, 1);
 
demodData_1          = (processedData0f_1 .^2 + processedData90f_1 .^2) .^(1/2);

%% Correct for amplitude of reference

demodDataC_1         = demodData_1*(2/modAmp_1);
     
meanF0               = mean(demodDataC_1);                                 % Mean value of baseline
medianF0             = median(demodDataC_1);                               % Median value of baseline
pc1F0                = prctile(demodDataC_1,1);                            % 1% percentile value of baseline 
pc5F0                = prctile(demodDataC_1,5);                            % 5% percentile value of baseline
pc10F0               = prctile(demodDataC_1,10);                           % 10% percentile value of baseline 
pc20F0               = prctile(demodDataC_1,20);                           % 20% percentile value of baseline
pc40F0               = prctile(demodDataC_1,40);                           % 40% percentile value of baseline 
pc80F0               = prctile(demodDataC_1,80);                           % 80% percentile value of baseline

%% F0 option
F0       = pc1F0;

%% Caclulate DFF
DFF                  = 100*((demodDataC_1-F0)/F0);

%% Get analog/digital AuxData from Oeps files - for comparison with behavior data
PhotometryData = [];
% adjust for clock offset between open ephys and kilosort
[offset] = AdjustClockOffset(myKsDir);
% adjust for clock offset between open ephys and kilosort
%timestamps = timestamps - offset; % for now timestamps is not right
Timestamps_LED_1 = Timestamps_LED_1 - offset; 

% downsample to behavior resolution
SampleRate = 500; % Samples/second

PhotometryData(:,1) = 0:1/SampleRate:max(Timestamps_LED_1); %timestamps based on behavior timestamps res
PhotometryData(:,2) = interp1q(Timestamps_LED_1,DFF,PhotometryData(:,1)); % DFF photometry signal with behavior timestamps res

% create a continuous TrialOn vector
for MyTrial = 1:size(TTLs.Trial,1)
    [~,start_idx] = min(abs(PhotometryData(:,1)-TTLs.Trial(MyTrial,1)));
    [~,stop_idx]  = min(abs(PhotometryData(:,1)-TTLs.Trial(MyTrial,2)));
    PhotometryData(start_idx:stop_idx,3) = 1;
end

%% create continuous trial on for odor 1???
for MyOdor = 1:size(TTLs.Odor1,1)
    [~,start_idx_od] = min(abs(PhotometryData(:,1)-TTLs.Odor1(MyOdor,1)));
    [~,stop_idx_od]  = min(abs(PhotometryData(:,1)-TTLs.Odor1(MyOdor,2)));
    PhotometryData(start_idx_od:stop_idx_od,4) = 1;
end

%% create vector with photometry data from 1s before trial start to 1s after
startoffset = 1;

for MyTrial = 1:size(TTLs.Trial,1)
    [~,start_idx] = min(abs(PhotometryData(:,1)-TTLs.Trial(MyTrial,1)));
    %[~,stop_idx]  = min(abs(PhotometryData(:,1)-TTLs.Trial(MyTrial,2)));
    begin_idx = start_idx - startoffset*SampleRate;
    end_idx = start_idx + startoffset*SampleRate;
    PhotometryTraces(:,MyTrial) = PhotometryData(begin_idx:end_idx,2);
end

%% Plots
figure()
subplot(3,1,1)
plot(PhotometryData(1:200000,2))
subplot(3,1,2)
plot(PhotometryData(1:200000,3))
subplot(3,1,3)
plot(PhotometryData(1:200000,4))

%%
figure()
plot(PhotometryData(1:200000,2)); hold on
plot(PhotometryData(1:200000,3)*3500)

%% 
figure()
PhotometryTrace_mean = mean(PhotometryTraces, 2);
PhotometryTrace_std = std(PhotometryTraces,0,2);
x = 1:1:size(PhotometryTrace_mean,1);
%shadedErrorBar(x,PhotometryTrace_mean,PhotometryTrace_std);
plot(PhotometryTrace_mean); hold on
xline(startoffset*SampleRate);

%% what session

switch whatmouse
    case 'DH2'
        WhichSession = 'DH2_20200312_r0_processed';
        SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\DH2';
        handles.WhereSession.String = fullfile(SessionPath,WhichSession);

    case 'DH4'
        WhichSession = 'DH4_20200310_r0_processed';
        SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\DH4';
        handles.WhereSession.String = fullfile(SessionPath,WhichSession);
end

%% Load relevant variables from behavior session

load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
handles.NumUnits.String = num2str(size(SingleUnits,2));
handles.SampleRate = SampleRate;


%% get the trial sorting order
allTrials = [];
for whichodor = 1:3
whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
    find(TrialInfo.Odor==whichodor));
OdorTrial = ones(size(whichTrials,1),1)*whichodor;
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) OdorTrial]; 
whichTrials = sortrows(whichTrials,2);
allTrials = [allTrials; whichTrials];
allTrials(find(allTrials(:,1)==1),:) = []; % getting rid of 1st trial 
end

%% 
for odor = 1:3
for TZ = 1:12
    whatTrials = allTrials(find(allTrials(:,2)==TZ & allTrials(:,3)==odor)); % find trials for this particular od and this particular TZ
    meanthisPhotometryTrace = mean(PhotometryTraces(:,whatTrials),2); % mean photometry signal for those specific trials
    stdthisPhotometryTrace = std(PhotometryTraces(:,whatTrials),0,2);
    meanPerTrialType(odor, TZ, :) = meanthisPhotometryTrace;
    stdPerTrialType(odor, TZ, :) = stdthisPhotometryTrace;
end    
end


%% Plots
figure()
x = 1:1:size(meanPerTrialType,3);
%i = 1;
for od = 1:3
    for tz = 1:12
        figure(od)
        subplot(4,3,tz)
        plot(squeeze(meanPerTrialType(od, tz, :)))
        %shadedErrorBar(x,meanPerTrialType(od, tz, :),stdPerTrialType(od, tz, :));
        xline(startoffset*SampleRate);
        %i = i+1;
    end
end

%% Get behavior traces from Trial x1 to x2
startTrial = 36; %x1
endTrial = startTrial +5; %x2

% load relevant variables
load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');
[TracesOut] = ConcatenateTraces(Traces, startTrial:endTrial, SampleRate*startoffset);

% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

% plot odor boxes on the behavior plot
% axes(handles.BehaviorPlot);
% hold off
figure()
subplot(2,1,1)
for i = 1:4
    handles.(['Trial',num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.(['Trial',num2str(i),'Plot']).EdgeColor = 'none';
    %ValveTS = TrialInfo.SessionTimestamps((TrialInfo.Odor(startTrial:endTrial)==i),1:2)' + TimestampAdjuster;  %changing to only include trials I want
    TS = TrialInfo.SessionTimestamps(startTrial:endTrial,1:2)+ TimestampAdjuster;
    ValveTS = TS((TrialInfo.Odor(startTrial:endTrial)==i),1:2)';
    %ValveTS = TrialInfo.SessionTimestamps((TrialInfo.Odor==i),1:2)' + TimestampAdjuster;
    if ~isempty(ValveTS)
        handles.(['Trial',num2str(i),'Plot']).Vertices = [ ...
            reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
            repmat([0 5 5 0]',size(ValveTS,2),1)];
        handles.(['Trial',num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
    end
    %set(gca,'XLim', [TracesOut.Timestamps{1, 1}(1) TracesOut.Timestamps{1, 1}(end)]);
end

% plot the target zone
handles.TargetZonePlot = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.2);
hold on;
handles.TargetZonePlot.EdgeColor = 'none';
TrialTS = TrialInfo.SessionTimestamps(startTrial:endTrial,1:2)' + TimestampAdjuster; %changing to only include trials I want
TZList =  TargetZones(TrialInfo.TargetZoneType(startTrial:endTrial),[3 1 1 3])'; %changing to only includ triald I want
handles.TargetZonePlot.Vertices = [ ...
        reshape([TrialTS(:) TrialTS(:)]', 2*numel(TrialTS), []) , ...
        TZList(:)];
handles.TargetZonePlot.Faces = reshape(1:2*numel(TrialTS),4,size(TrialTS,2))';
    
% plot the lever trace on top
plot(Timestamps + TimestampAdjuster, TracesOut.Lever{1},'k');

% plot Rewards
handles.reward_plot = plot(NaN, NaN, 'color',Plot_Colors('t'),'Linewidth',1.25);
tick_timestamps =  Timestamps(find(diff(TracesOut.Rewards{1}==1)) + 1)' + TimestampAdjuster;
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; 5; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.reward_plot,'XData',tick_x,'YData',tick_y);


set(gca,'XLim', [TracesOut.Timestamps{1, 1}(1)+TimestampAdjuster TrialInfo.SessionTimestamps(endTrial,2)+1+TimestampAdjuster]);
%end

%% Extract Trials from Trial x1 to x2

[~,start_idx] = min(abs(PhotometryData(:,1)-TTLs.Trial(startTrial,1))); % timestamps of x1 trial in behavior stamps
[~,stop_idx]  = min(abs(PhotometryData(:,1)-TTLs.Trial(endTrial,2))); % timestamps of x2 trial in behavior stamps
begin_idx = start_idx - startoffset*SampleRate; %-1s
end_idx = stop_idx + startoffset*SampleRate; %%+1s
myPhotometryTime = PhotometryData(begin_idx:end_idx,1);
myPhotometryTrace = PhotometryData(begin_idx:end_idx,2);

subplot(2,1,2)
for i = 1:4
    handles.(['Trial',num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.(['Trial',num2str(i),'Plot']).EdgeColor = 'none';
    TS = TrialInfo.SessionTimestamps(startTrial:endTrial,1:2)+ TimestampAdjuster;
    ValveTS = TS((TrialInfo.Odor(startTrial:endTrial)==i),1:2)';
    if ~isempty(ValveTS)
        handles.(['Trial',num2str(i),'Plot']).Vertices = [ ...
            reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
            repmat((10+max(myPhotometryTrace))*[0 1 1 0]',size(ValveTS,2),1)];
        handles.(['Trial',num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
    end
    %set(gca,'XLim', [TracesOut.Timestamps{1, 1}(1) TracesOut.Timestamps{1, 1}(end)]);
end
plot(myPhotometryTime,myPhotometryTrace, 'k')
set(gca, 'XLim', [TracesOut.Timestamps{1, 1}(1)+TimestampAdjuster TrialInfo.SessionTimestamps(endTrial,2)+1+TimestampAdjuster],'YLim', [(min(myPhotometryTrace)-10) (10+max(myPhotometryTrace))])



%% WIP
% TracePerTrialType = [];
% for odor = 1:3
% for TZ = 1:12
%     whatTrials = allTrials(find(allTrials(:,2)==TZ & allTrials(:,3)==odor));
%     thisPhotometryTrace = PhotometryTraces(:,whatTrials);
%     TracePerTrialType(odor, TZ, :) = [TracePerTrialType thisPhotometryTrace];
% end    
% end