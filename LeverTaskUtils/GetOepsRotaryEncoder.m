function [Position, FrameTriggers] = GetOepsRotaryEncoder(myKsDir, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('SampleRate', 500, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
SampleRate = params.Results.SampleRate;

%% Get Encoder clicks from the OpenEphys Events file
if exist(fullfile(myKsDir,'myTTLfile_1.mat'))
    foo             = load(fullfile(myKsDir,'myTTLfile_1.mat'),'TTLs');
    data            = foo.TTLs.data;
    timestamps      = foo.TTLs.timestamps;
    info.eventId    = foo.TTLs.info.eventId;
    offset          = foo.TTLs.offset;
else
    FrameTriggers = [];
    Position = [];
    disp('no TTLs file found in the sorting folder');
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

% create an analog position vector with continuous timestamps
% downsampled to behavior resolution
Position(:,1) = 0:1/SampleRate:EncoderEvents(end,1);
Position(:,2) = interp1q(EncoderEvents(:,1),EncoderEvents(:,4),Position(:,1));

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
if ~isempty(timestamps(intersect(find(~info.eventId),find(data==8))))
    VideoTS(:,2) = timestamps(intersect(find(~info.eventId),find(data==8)));
else
    VideoTS(end,2) = timestamps(end);
end

for i = 1:size(VideoTS,1)
    whichFrames = find(FrameTriggers(:,1)>VideoTS(i,1) & FrameTriggers(:,2)<VideoTS(i,2));
    FrameTriggers(whichFrames,4) = i;
end

end
