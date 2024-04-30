SessionName = 'Q4_20221109_r0';
MySession = fullfile('/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/', ...
                        SessionName(1:regexp(SessionName,'_','once')-1), ...
                        [SessionName,'_processed.mat']);
load(MySession,'Traces','TimestampAdjust');
TracesOut = ConcatenateTraces2Mat(Traces);

% plot the respiration trace (in OEPS base)
plot(TracesOut(:,8) + TimestampAdjust.ClosedLoop ,TracesOut(:,1))

fbinary = '/mnt/albeanu_lab/priyanka/EphysData/Q4/2022-11-09_13-04-41';
% read the open ephys flat binary file
session = Session(fbinary);
TotalSamples = size(session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps,1);

TS = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps;
TS = TS - session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps(1);

Lever = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(43,:);

% do lever from behavior and OEPS match
plot(TS,2.5 + (double(Lever)/12800));
hold on
plot(TracesOut(:,8) + TimestampAdjust.ClosedLoop ,TracesOut(:,1));

% get a single unit channel : ch33
Electrode = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(33,:)';

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

