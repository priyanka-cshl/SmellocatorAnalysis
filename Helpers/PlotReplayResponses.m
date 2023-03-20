function [] = PlotReplayResponses(SessionPath,MyUnits)

if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
elseif strcmp(computer,'PCWIN64')
    datapath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
elseif strcmp(computer, 'GLNXA64') % nicobar
    datapath = '/mnt/data/Processed/Behavior/';
else
    datapath = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

%% which Units to use
N = size(SingleUnits,2); % #units
if nargin < 2 || isempty(MyUnits)
    MyUnits = 1:N;
    % sort units by tetrodes
    foo = cell2mat(arrayfun(@(x) [x.tetrode; x.id], SingleUnits, 'UniformOutput', false))';
    [~, MyUnits] = sortrows(foo,[1 2]);
end


%% Get the replay traces and spikes
[MyTraces,timestamps,PSTH,Raster] = ...
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
                        'plotfigures', 1, 'plotephys', 1, ...
                        'PlotOpenLoop', 1, ...
                            'whichunits', MyUnits, 'UnitsPerFig', 5);

end % function end


