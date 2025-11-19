function [TTLs,AuxData] = GetOepsRotaryEncoder(myKsDir, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('Analog', true, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
GetPosition = params.Results.Analog;

%% defaults
SampleRate = 500; % Behavior Acquisition rate - standard rate at which we downsample analog data
OEPSSamplingRate = 30000; % needed to read the aux file
PPR = 1024; % pulses per rotation for the encoder
angularStep = (1024*4)^-1; % by 4 for the quadrature encoder

%% Get Encoder clicks from the OpenEphys Events file
if exist(fullfile(myKsDir,'myTTLfile_1.mat'))
    foo             = load(fullfile(myKsDir,'myTTLfile_1.mat'),'TTLs');
    data            = foo.TTLs.data;
    timestamps      = foo.TTLs.timestamps;
    info.eventId    = foo.TTLs.info.eventId;
    offset          = foo.TTLs.offset;
else
    TTLs = [];
    AuxData = [];
    disp('no TTLs found in the sorting folder');
    return
end

timestamps = timestamps - offset;

%% Get all Encoder Events
whichEvents = find(data==3 | data==5);
EncoderEvents = [timestamps(whichEvents) double(info.eventId(whichEvents)) double(data(whichEvents))];
EncoderEvents = sortrows(EncoderEvents,1); % sort by time

EncoderEvents(:,4) = nan;
EncoderEvents(1,4) = 1;
for i = 2:size(EncoderEvents,1)
    if EncoderEvents(i,3) ~= EncoderEvents(i-1,3)
        EncoderEvents(i,4) = EncoderEvents(i-1,4) + 1;
    else
        EncoderEvents(i,4) = EncoderEvents(i-1,4) - 1;
    end

end

%% Get Camera triggers
On = timestamps(intersect(find(info.eventId),find(data==6)));
Off = timestamps(intersect(find(~info.eventId),find(data==6)));
if Off(1)<On(1) % first off value, preceeds the first On
    On = vertcat(nan, On);
end
if On(end)>Off(end) % last on value, exceeds the last Off
    Off = vertcat(Off, nan);
end
FrameTriggers = [On Off Off-On];
% find out which triggers were in the recorded video
VideoTS(:,1) = timestamps(intersect(find(info.eventId),find(data==8)));
VideoTS(:,2) = timestamps(intersect(find(~info.eventId),find(data==8)));

for i = 1:size(VideoTS,1)
    whichFrames = find(FrameTriggers(:,1)>VideoTS(i,1) & FrameTriggers(:,2)<VideoTS(i,2));
    FrameTriggers(whichFrames,4) = i;
end

%% Read the AuxFile
fid = fopen(fullfile(myKsDir,'myauxfile.dat'));
AuxData = fread(fid,'*int16');

% get no. of samples from session details
load(fullfile(myKsDir,'SessionDetails.mat'),'Files');
NSamples    = Files.Samples;
nAuxChans   = numel(AuxData)/NSamples;
plot(AuxData(4:8:end));


OEPSSamplingRate = handles.OEPSSamplingRate; %30000;
Nchan       = str2double(handles.NumChans.String);
Nsamples    = str2double(handles.WindowSize.String)*OEPSSamplingRate;

TStart      = handles.TimeScroll.Value * 60 * ...
                str2double(handles.RecordingLength.String) - ...
                str2double(handles.WindowSize.String);
TStart = floor(TStart);
offset = uint64(max(0,2*Nchan * TStart * OEPSSamplingRate)); % 2 bytes per sample

fid = fopen(fullfile(handles.BinaryPath.String,'mybinaryfile.dat'),'r');
fseek(fid, offset, 'bof');
MyData = fread(fid, [Nchan Nsamples], '*int16');

AuxData = [];
if GetAux
    if UseSortingFolder
        %if exist(fullfile(myKsDir,'myTTLfile_1.mat'))
        
    else

        %% Get analog/digital AuxData from Oeps files - for comparison with behavior data
        foo = dir(fullfile(myKsDir,'*_ADC1.continuous')); % pressure sensor
        filename = fullfile(myKsDir,foo.name);
        [Auxdata1, timestamps, ~] = load_open_ephys_data(filename); % data has channel IDs
        foo = dir(fullfile(myKsDir,'*_ADC2.continuous')); % thermistor
        filename = fullfile(myKsDir,foo.name);
        [Auxdata2, ~, ~] = load_open_ephys_data(filename); % data has channel IDs

        % adjust for clock offset between open ephys and kilosort
        timestamps = timestamps - offset;

        % downsample to behavior resolution

        AuxData(:,1) = 0:1/SampleRate:max(timestamps);
        AuxData(:,2) = interp1q(timestamps,Auxdata1,AuxData(:,1)); % pressure sensor
        AuxData(:,3) = interp1q(timestamps,Auxdata2,AuxData(:,1)); % thermistor
        % create a continuous TrialOn vector
        for MyTrial = 1:size(TTLs.Trial,1)
            [~,start_idx] = min(abs(AuxData(:,1)-TTLs.Trial(MyTrial,1)));
            [~,stop_idx]  = min(abs(AuxData(:,1)-TTLs.Trial(MyTrial,2)));
            AuxData(start_idx:stop_idx,4) = 1;
        end
    end
end
end
