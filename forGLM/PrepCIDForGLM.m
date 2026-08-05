function [TracesOut,SingleUnits] = PrepCIDForGLM(myKsDir,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('whichSniffSensor', 'mfs', @(x) ischar(x)); % 'mfs', 'therm', 'both'
params.addParameter('allUnits', 0, @(x) isnumeric(x)); % -1: non-good, 0: good (default), 1: all
params.addParameter('minRate', 0, @(x) isnumeric(x)); % in Hz
params.addParameter('savemode', 0, @(x) isnumeric(x));
params.addParameter('sampleRate', 500, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
whichSniffSensor = 1+ strcmp(params.Results.whichSniffSensor,'mfs'); % 1 = thermistor, 2 = MFS, 3 = both
UnitFilter = params.Results.allUnits; % -1: non-good, 0: good (default), 1: all
minRate = params.Results.minRate; % in Hz
savemode = params.Results.savemode;
sampleRate = params.Results.sampleRate; % in Hz

%% Load the relevant data
[StimSettings, TTLs, SingleUnits, ~, TrialWiseSniffs, ~] = ...
    CIDResponsePrepper(myKsDir,'whichSniffSensor',params.Results.whichSniffSensor,'allUnits',UnitFilter,'minRate',minRate);
% PID kernels
PIDfit = load('/mnt/data/CID/PID/16odor_basic_fit.mat');
PIDdt = mode(diff(PIDfit.F.TimeVector));

% Sniffs
SniffDeltas = TrialWiseSniffs(:,1) + TrialWiseSniffs(:,6); % inhalation starts
% the first sniff is missing here - add that
SniffDeltas = [(SniffDeltas(1) - TrialWiseSniffs(1,8)); SniffDeltas];

% Make a time trace
maxT = SniffDeltas(end) + 1; % 1 second after the last sniff
Timestamps = (0:(1/sampleRate):maxT)';
ZeroBasis = Timestamps*0;

% make a digital sniff trace
DigitalSniffs = ZeroBasis;
for i = 1:numel(SniffDeltas)
    [~,idx] = min(abs(Timestamps-SniffDeltas(i)));
    DigitalSniffs(idx) = 1;
end

if strcmp(StimSettings.SessionType, '16_Odors')
    % TTLs to odor waveforms
    DigitalOdor = ZeroBasis;
    PIDOdor     = ZeroBasis;
    for i = 1:size(TTLs.Trial,1) % every trial
        % only use the first odor pulse
        TS          = TTLs.Trial(i,[7 8]); % this trial odor ON and OFF
        if i == 1
            startAt = find(Timestamps>=TTLs.Trial(i,7),1,'first');
            endAt   = find(Timestamps<TTLs.Trial(i+1,7),1,'last'); % until the next odor pulse
        elseif i == size(TTLs.Trial,1)
            endAt = size(Timestamps,1);
        else
            endAt   = find(Timestamps<TTLs.Trial(i+1,7),1,'last'); % until the next odor pulse
        end
        valveState = [];
        valveState(:,1) = Timestamps(startAt:endAt); % from odor On to next odor ON
        odorIdx = find((valveState>=TS(1))&(valveState<=TS(2)));
        valveState(:,2) = valveState(:,1)*0;
        valveState(odorIdx,2) = 1;
        whichOdor = TTLs.Trial(i,4);
        PIDtrace = ValveToPIDOdor(valveState,PIDfit.F.paramsout(:,whichOdor),PIDdt);
        % normalize
        PIDtrace = PIDtrace/max(PIDtrace);
        PIDOdor(startAt:endAt) = PIDtrace;

        DigitalOdor(startAt:endAt) = -whichOdor;
        % map indices to actual timestamps
        DigitalOdor((startAt-1)+odorIdx) = whichOdor;

        startAt = endAt + 1;
    end
end

%%
TracesOut.Timestamps{1} = Timestamps;
TracesOut.SniffsDigitized{1} = DigitalSniffs;
TracesOut.Odor{1} = DigitalOdor;
TracesOut.PID{1} = PIDOdor;

%% Save to Disk
if savemode
    [~,MouseName] = fileparts(fileparts(myKsDir));
    [~,filename] = fileparts(myKsDir);
    filename = [MouseName,'_',regexprep(filename(1,1:10),'-',''),'_r0_processed.mat'];
    savepath = '/mnt/data/';
    save(fullfile(savepath,'forWDW',filename),'TracesOut','SingleUnits');
end

end
