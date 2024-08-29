% script to compare sniff peak/valley detections in thermistor vs. mass
% flow sensor

Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));

DataRoot = '/mnt/grid-hs/mdussauz/SniffingTest/2023-07-07_17-01-58_OEPS_format/Record Node 102';
MFfile = '100_RhythmData_ADC1.continuous';
THfile = '100_RhythmData_ADC3.continuous';

% load both the sensors, downsample and filter
[MFraw, Auxtimestamps, ~] = load_open_ephys_data(fullfile(DataRoot,MFfile));
[THraw, Auxtimestamps, ~] = load_open_ephys_data(fullfile(DataRoot,THfile));

Auxtimestamps = (Auxtimestamps - Auxtimestamps(1))/30000;

% downsample to 500hz
SniffData = [];
SniffData(:,1) = 0:1/500:max(Auxtimestamps);
SniffData(:,2) = interp1q(Auxtimestamps,MFraw,SniffData(:,1));
SniffData(:,3) = FilterThermistor(SniffData(:,2));
SniffData(:,4) = interp1q(Auxtimestamps,THraw,SniffData(:,1));
SniffData(:,5) = FilterThermistor(SniffData(:,4));

% detect peaks/valleys in the thermistor data
[TH_TS,TH_Indices] = ProcessThermistorData(SniffData(:,[1 5]));