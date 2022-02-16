function [Corrs,g,Residuals,Tags] = CompareReplayPSTHs(SessionPath,MyUnits,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 50, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
Binsize = params.Results.binsize;

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

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

BehaviorBinsize = Binsize/(1000/SampleRate);

%% which Units to use
N = size(SingleUnits,2); % #units
if nargin < 2 || isempty(MyUnits)
    MyUnits = 1:N;
end

%% Get the replay traces and spikes
[OpenLoopTraces,OpenLoopTimestamps,OpenLoopPSTH,OpenLoopRasters] = ...
    ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, 'whichunits', MyUnits);

Trialcounts = [1 numel(find(ReplayTTLs.TrialID<=TrialInfo.TrialID(end))) numel(find(ReplayTTLs.TrialID>TrialInfo.TrialID(end)))];

%% Parse the continuous traces into Trials
Trial = OpenLoopTraces(:,6); % Trial vector
Trial(Trial<0) = 0; % ignore OdorOn periods
Odor = abs(OpenLoopTraces(:,6)); % Trial + OdorON
TrialTS =  horzcat( find(diff(Odor)>0), find(diff(Trial)>0), find(diff(Trial)<0)); % Odor ON, Trial ON, Trial OFF
TrialTS(:,4) = Odor(TrialTS(:,2)); % which odor
TrialTS(:,5) = OpenLoopTraces(TrialTS(:,2),7); % which Target Zone


for t = 1:size(TrialTS,1) % every subtrial
    % get which trace indices to use
    idx(1) = TrialTS(t,1); % odor start
    if t < size(TrialTS,1)
        idx(2) = TrialTS(t+1,1) - 1; % next trial odor start
    else
        idx(2) = TrialTS(t,3) + SampleRate; % 1 sec post trial off
    end
    % adjust idx(2) to allow for binning
    idx(2) = idx(1) + BehaviorBinsize*floor(numel(idx(1):idx(2))/BehaviorBinsize) - 1;
    idx = idx * (Binsize/BehaviorBinsize); % Spike rasters are in ms resolution
    idx(1) = idx(1) - (Binsize/BehaviorBinsize) + 1;
    BinnedPSTH = [];
    for i = 1:numel(MyUnits)
        for rep = 1:size(OpenLoopRasters,1) % every repeat
            thisRepPSTH = sum(reshape(OpenLoopRasters(rep,idx(1):idx(2),i),Binsize,[]));
            BinnedPSTH(rep,:,i) = thisRepPSTH;
        end
    end
    [Corrs(:,:,t),g] = ReplayCrossCorr(BinnedPSTH,Trialcounts);
    [Residuals(:,:,t),Tags] = ReplayResiduals(BinnedPSTH,Trialcounts);
end


end % function end


