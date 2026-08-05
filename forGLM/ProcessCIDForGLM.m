function [TracesOut,SingleUnits] = ProcessCIDForGLM(myKsDir,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('whichSniffSensor', 'mfs', @(x) ischar(x)); % 'mfs', 'therm', 'both'
params.addParameter('allUnits', 0, @(x) isnumeric(x)); % -1: non-good, 0: good (default), 1: all
params.addParameter('minRate', 0, @(x) isnumeric(x)); % in Hz
params.addParameter('savemode', 0, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
whichSniffSensor = 1+ strcmp(params.Results.whichSniffSensor,'mfs'); % 1 = thermistor, 2 = MFS, 3 = both
UnitFilter = params.Results.allUnits; % -1: non-good, 0: good (default), 1: all
minRate = params.Results.minRate; % in Hz
savemode = params.Results.savemode;

%WhereSession = '/mnt/data/Sorted/T3/2025-05-08_16-20-35/quickprocesssniffs.mat';
WhereSession = fullfile(myKsDir,'quickprocesssniffs.mat');

%% Load the relevant data
[StimSettings, TTLs, SingleUnits, ~, TrialWiseSniffs, ~] = ...
    CIDResponsePrepper(myKsDir,'whichSniffSensor',params.Results.whichSniffSensor,'allUnits',UnitFilter,'minRate',minRate);

% Sniffs
SniffDeltas = TrialWiseSniffs(:,1) + TrialWiseSniffs(:,6); % inhalation starts
% the first sniff is missing here - add that
SniffDeltas = [(SniffDeltas(1) - TrialWiseSniffs(1,8)); SniffDeltas];

% Spikes are in SingleUnits.spikes

% PID kernels
PIDfit = load('/mnt/data/CID/PID/16odor_basic_fit.mat');


% make a time-series with sniff onsets
DigitalSniffs = TracesOut.SniffsFiltered{1}*0;
if size(CuratedSniffTimestamps,2) < 10
    CuratedSniffTimestamps(:,10) = 0;
end
for n = 1:size(CuratedSniffTimestamps,1)
    idx = CuratedSniffTimestamps(n,8:9);
    if CuratedSniffTimestamps(n,10) == -1
        DigitalSniffs(idx(1):idx(2)) = -1;
    else
        DigitalSniffs(idx(1):idx(2)) = 1;
    end
end
TracesOut.SniffsDigitized{1} = DigitalSniffs;







% Sniffs
[SniffCoords,SniffProps] = TallyThermistorNPressureSniffs(WhereSession);
if ~isempty(SniffCoords)
    load(WhereSession, 'RespirationData');
else
    try
        load(WhereSession, 'CuratedMFSSniffTimestamps','RespirationData');
        CuratedSniffTimestamps = CuratedMFSSniffTimestamps;
    catch
        load(WhereSession, 'CuratedSniffTimestamps','RespirationData');
    end
end

try
    load(fileparts(WhereSession),'quickprocessOdorTTLs.mat');
    if ~exist('TTLs','var')
        TTLs = [];
    end
catch
    TTLs = [];
end

%% Load units
SingleUnits = LoadKS4Units(myKsDir,'minSpikes',minRate,'allUnits',UnitFilter);
nUnits = size(SingleUnits,2);

%% Traces
TracesOut.Timestamps{1}         = RespirationData(:,1);
% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';
% snity check for timestamp drops
if any(round(diff(Timestamps),3,'decimal')~=0.002)
    keyboard;
end

if isempty(TTLs)
    TracesOut.Odor{1} = TracesOut.Timestamps{1}*0;
    TracesOut.Trial{1} = TracesOut.Timestamps{1}*0;
    TracesOut.Manifold{1} = TracesOut.Timestamps{1}*0;
else
    % stimulus identity
    % set sniffs we don't want to use to -1 or something like that
    % air sniffs = TracesOut.Manifold{1} = 1
%     TS          = [0 TTLs.Trial(1,1)];
%     Samps       = intersect(find(TracesOut.Timestamps{1}>=TS(1)),find(TracesOut.Timestamps{1}<=TS(2)));
%     TracesOut.Manifold{1}(Samps) = 1; 
    TracesOut.Odor{1} = TracesOut.Timestamps{1}*0;
    TracesOut.Trial{1} = TracesOut.Timestamps{1}*0;
    TracesOut.Manifold{1} = 1 + TracesOut.Timestamps{1}*0; % air is ON all throughout in CID experiments

    % odor sniffs = TracesOut.Odor{1} == stimID

    for i = 1:size(TTLs.Trial,1) % every trial
        % only use the first odor pulse
        TS          = TTLs.Trial(i,[1 2 7 8]);
        Samps       = intersect(find(TracesOut.Timestamps{1}>=TS(1)),find(TracesOut.Timestamps{1}<=TS(2)));
        TracesOut.Trial{1}(Samps) = 1;
        Samps       = intersect(find(TracesOut.Timestamps{1}>=TS(3)),find(TracesOut.Timestamps{1}<=TS(4)));
        TracesOut.Odor{1}(Samps) = TTLs.Trial(i,4);
    end
end

% fake trace for fitting kernels similar to that for the lever task
TracesOut.Motor{1} = TracesOut.Timestamps{1}*0;

%% sniffing specific
% add a filtered sniff trace
if any(RespirationData(:,3))
    TracesOut.SniffsFiltered{1}     = RespirationData(:,3);
    if size(RespirationData,2) == 4
        TracesOut.MassFlowSensor{1}     = RespirationData(:,4) - 2.5;
    end
elseif any(RespirationData(:,4))
    MFS      = RespirationData(:,4) - 2.5; % unfiltered MassFlowSensor trace
    % filter
    fband = [0.1 30];
    Np    = 4; % filter order
    [b,a] = butter(Np,fband/(500/2));
    % predicted thermistor trace
    MFS2Therm= filtfilt(b,a,cumsum(MFS));
    TracesOut.SniffsFiltered{1} = MFS2Therm;
    TracesOut.MassFlowSensor{1} = MFS;
else
    keyboard;
end

% fake trace for fitting kernels similar to that for the lever task
TracesOut.SniffsLocationed{1}   = TracesOut.SniffsFiltered{1}*nan;

if ~isempty(SniffCoords)
    % make several digitized sniffs
    DigitalSniffs = zeros(numel(TracesOut.SniffsFiltered{1}),3);
    for n = 1:size(SniffCoords,1)
        idx = SniffCoords(n,[4 5 9 10 14 15]);
        idx = reshape(idx,2,[])';
        for m = 1:3
            if ~any(isnan(idx(m,:)))
                if SniffCoords(n,1)>0
                    DigitalSniffs(idx(m,1):idx(m,2),m) = 1;
                else
                    DigitalSniffs(idx(m,1):idx(m,2),m) = -1;
                end
            end
        end
    end
    TracesOut.SniffsDigitized{1} = DigitalSniffs(:,1); % thermistor peaks
    TracesOut.MFS2ThermDigitized{1} = DigitalSniffs(:,2); % mfs to thermistor peaks
    TracesOut.MFSDigitized{1} = DigitalSniffs(:,3); % mfs zero crossings
else
    DigitalSniffs = TracesOut.SniffsFiltered{1}*0;
    if size(CuratedSniffTimestamps,2) < 10
        CuratedSniffTimestamps(:,10) = 0;
    end
    for n = 1:size(CuratedSniffTimestamps,1)
        idx = CuratedSniffTimestamps(n,8:9);
        if CuratedSniffTimestamps(n,10) == -1
            DigitalSniffs(idx(1):idx(2)) = -1;
        else
            DigitalSniffs(idx(1):idx(2)) = 1;
        end
    end
    TracesOut.SniffsDigitized{1} = DigitalSniffs;
end

%% wheel position data and camera triggers if video was recorded
[WheelPosition, FrameTriggers] = GetOepsRotaryEncoder(myKsDir);
if ~isempty(WheelPosition)
    if size(WheelPosition,1) <=  size(TracesOut.Timestamps{1},1)
        nWheel = size(WheelPosition,1);
        nanWheel = find(~isnan(WheelPosition(:,2)),1,"first");
        WheelPosition(1:nanWheel,2) = WheelPosition(nanWheel,2);
        WheelPosition((nWheel+1):size(TracesOut.Timestamps{1},1),2) = WheelPosition(nWheel,2);
        TracesOut.WheelPosition{1} = WheelPosition(:,2);
    else
        keyboard;
    end

    if any(FrameTriggers(:,4)) % video was recorded
        % ignore nans
        nanRows = find(isnan(FrameTriggers(:,1))|isnan(FrameTriggers(:,2)));
        FrameTriggers(nanRows,:) = [];
        fakeTriggers = find(FrameTriggers(:,3)<0.2*mode(FrameTriggers(:,3)));
        FrameTriggers(fakeTriggers,:) = [];
        TracesOut.VideoTriggers{1} = TracesOut.Timestamps{1}*0;
        for i = 1:size(FrameTriggers,1)
            [~,idx1] = min(abs(TracesOut.Timestamps{1}-FrameTriggers(i,1)));
            [~,idx2] = min(abs(TracesOut.Timestamps{1}-FrameTriggers(i,2)));
            if FrameTriggers(i,4)
                TracesOut.VideoTriggers{1}(idx1:idx2) = 1;
            else
                TracesOut.VideoTriggers{1}(idx1:idx2) = -1;
            end
            %FrameTriggers(i,5:6) = [idx1 idx2];
        end
    end
end

%% Save to Disk
if savemode
    [~,MouseName] = fileparts(fileparts(myKsDir));
    [~,filename] = fileparts(myKsDir);
    filename = [MouseName,'_',regexprep(filename(1,1:10),'-',''),'_r0_processed.mat'];
    savepath = '/mnt/data/';
    save(fullfile(savepath,'forWDW',filename),'TracesOut','SingleUnits');
end

end
