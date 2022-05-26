function [OpenLoop] = ProcessReplayResponses(MySession,MyUnits)

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

%% which Units to use
N = size(SingleUnits,2); % #units
if nargin < 2 || isempty(MyUnits)
    MyUnits = 1:N;
end

%% Get the replay traces and spikes
[MyTraces,timestamps,PSTH,Raster] = ...
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
                        'plotfigures', 0, 'plotephys', 0, ...
                        'PlotOpenLoop', 0, ...
                            'whichunits', MyUnits);
                        
                        
%% bunch of stuff to simplify things for giving data to Alex

% make all reward vectors identical
temprewards = sum(MyTraces(:,5,:),3);
temprewards(temprewards<2) = 0;
temprewards(temprewards>0) = 1;
for i = 1:size(MyTraces,3)
    MyTraces(:,5,i) = temprewards;
end

clear OpenLoop

OpenLoop.Traces.Lever   = squeeze(MyTraces(:,1,:));
OpenLoop.Traces.Motor   = squeeze(MyTraces(:,2,:));
OpenLoop.Traces.Sniffs  = squeeze(MyTraces(:,3,:));
OpenLoop.Traces.Rewards = squeeze(MyTraces(:,5,:));
OpenLoop.Traces.Trial   = squeeze(MyTraces(:,6,:));
OpenLoop.Traces.Timestamps = timestamps;


% chop PSTH and Rasters longer than the behavior data
idx = size(MyTraces,1)+1;
if size(PSTH,2)>idx
    PSTH(:,idx:end,:) = [];
    Raster(:,(1+(idx-1)*2):end,:) = [];
end
OpenLoop.PSTH = PSTH;
OpenLoop.Rasters = Raster;

OpenLoop.Trialcounts = [1 size(MyTraces,3)-1 size(Raster,1)];
OpenLoop.Trialcounts(3) = OpenLoop.Trialcounts(3) - 1 - OpenLoop.Trialcounts(2);

end % function end


