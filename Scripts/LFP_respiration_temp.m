SessionName = 'Q4_20221109_r0';
SessionName = 'O3_20211005_r0';
SessionName = 'Q8_20221204_r0';
MySession = fullfile('/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/', ...
                        SessionName(1:regexp(SessionName,'_','once')-1), ...
                        [SessionName,'_processed.mat']);
load(MySession,'Traces','TimestampAdjust');
TracesOut = ConcatenateTraces2Mat(Traces);

% plot the respiration trace (in OEPS base)
plot(TracesOut(:,7) + TimestampAdjust.ClosedLoop ,TracesOut(:,1))

fbinary = '/mnt/albeanu_lab/priyanka/EphysData/Q4/2022-11-09_13-04-41';
fbinary = '/mnt/grid-hs/mdussauz/ephysdata/lever_task/BatchO/O3/2021-10-05_14-24-31/Record Node 108';
fbinary = '/mnt/albeanu_lab/priyanka/EphysData/Q8/2022-12-04_22-40-13';
fbinary = '/mnt/albeanu_lab/priyanka/EphysData/O9/2022-06-30_15-14-32';

%% read the open ephys flat binary file
session = Session(fbinary);
TotalSamples = size(session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps,1);

TS = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps;
TS = TS - session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps(1);

%%
% filter in LFP band
fshigh = 0.1; fslow = 300; fs = 30000;
[b, a] = butter(3, [fshigh/fs,fslow/fs]*2, 'bandpass');

HighCutOff = 500;
LowCutOff = 20;
[b_high, a_high] = butter(5,HighCutOff/fs,"high"); % filter co-efficients
[b_low, a_low] = butter(5,LowCutOff/fs,"low");
d = designfilt('bandstopiir', 'FilterOrder',2,...
    'HalfPowerFrequency1',40,'HalfPowerFrequency2',70,...
    'DesignMethod','butter','SampleRate',fs);
myxlim = [50 60];

figure; subplot(3,1,1); hold on; subplot(3,1,2); hold on
indices = 30000 + (1:(30000*60));
for i = 1:4 %0
    neural.raw = double(session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(i,indices))';
    % filter the data
    neural.spikes = filtfilt(b_high,a_high,neural.raw); % filter the signals using zero phase filtering
    neural.LFP = filtfilt(b_low,a_low,neural.raw);
    neural.notched = filtfilt(d,neural.LFP);
        
%     ch_filtered = filter(b,a,raw);
%     ch_filtered = flipud(ch_filtered);
%     ch_filtered = filter(b,a,ch_filtered);
%     ch_filtered = flipud(ch_filtered);
%     subplot(3,1,1);
%     if mod(i,2)
%         plot(TS(indices),(i-1)+(neural.spikes/1000),'r');
%     else
%         plot(TS(indices),(i-1)+(neural.spikes/1000),'k');
%     end
    
    subplot(3,1,1);
    if mod(i,2)
        plot(TS(indices),(i-1)+(neural.notched/1000),'r');
    else
        plot(TS(indices),(i-1)+(neural.notched/1000),'k');
    end
    
    subplot(3,1,2);
    if mod(i,2)
        plot(TS(indices),(i-1)+(neural.LFP/1000),'r');
    else
        plot(TS(indices),(i-1)+(neural.LFP/1000),'k');
    end
end
set(gca,'XLim',myxlim);
subplot(3,1,1)
plot(TracesOut(:,7) + TimestampAdjust.ClosedLoop ,(TracesOut(:,3)*6)+5);
set(gca,'XLim',myxlim);

% plot thermistor
subplot(3,1,3);
plot(TracesOut(:,7) + TimestampAdjust.ClosedLoop ,(TracesOut(:,3)*6)+5);
set(gca,'XLim',myxlim);
%%

Lever = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(43,:);

% do lever from behavior and OEPS match
plot(TS,2.5 + (double(Lever)/12800));
hold on
plot(TracesOut(:,8) + TimestampAdjust.ClosedLoop ,TracesOut(:,1));

% get a single unit channel : ch33
Electrode = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(33,:)';
fs = 30000;
Chunk = double(Electrode(1:(fs*60))); % one minute of data
dt = 1/fs;
plot(dt:dt:60,Chunk)
subplot(2,1,1); plot(dt:dt:60,Chunk)

% low-pass filter
fcut = 300;
[b, a] = butter(3, [fcut/fs]*2, 'low');

% filter the data
Filtered = filter(b,a,Chunk);
Filtered = flipud(Filtered);
Filtered = filter(b,a,Filtered);
Filtered = flipud(Filtered);

% bandstop the AC band
Flt = bandstop(Filtered,[48 52], fs);


lpFilt = designfilt('lowpassiir','FilterOrder',16, ...
    'PassbandFrequency',10,'PassbandRipple',0.2, ...
    'SampleRate',30000);


fs = 30000;
d = designfilt('bandpassiir','FilterOrder',2, ...
               'HalfPowerFrequency1',2,'HalfPowerFrequency2',30, ...
               'DesignMethod','butter','SampleRate',fs);
           
           lpFilt = designfilt('lowpassiir','FilterOrder',16, ...
         'PassbandFrequency',10,'PassbandRipple',0.2, ...
         'SampleRate',30000);

E_f = filtfilt(d,Electrode);


Ef = bandstop(Electrode,[59 70],fs);

% filter in standard ephys band
fshigh = 300; fslow = 6000; fs = 30000;
[b, a] = butter(3, [fshigh/fs,fslow/fs]*2, 'bandpass');

% filter the data
Electrode = filter(b,a,Electrode);
Electrode = flipud(Electrode);
Electrode = filter(b,a,Electrode);
Electrode = flipud(Electrode);

% get a single unit channel : ch33
Electrode = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(33,:)';

% filter in LFP band
fshigh = 0.1; fslow = 30; fs = 30000;
[b, a] = butter(3, fslow/fs*2, 'low');

fshigh = 0.1; fslow = 30; fs = 30000;
[b, a] = butter(3, [fshigh/fs,fslow/fs]*2, 'bandpass');

% filter the data
LFP = filter(b,a,Electrode);
LFP = flipud(LFP);
LFP = filter(b,a,LFP);
LFP = flipud(LFP);

% plot thermistor (from behavior) and LFP
plot(TS,LFP/30);
hold on
plot(TracesOut(:,8) + TimestampAdjust.ClosedLoop ,TracesOut(:,3));


% Thermistor = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(42,:);

