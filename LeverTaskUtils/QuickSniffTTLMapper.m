%function [AllSniffs] = QuickSniffTTLMapper(myKsDir)

%myKsDir = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Ephys/Q9/2022-11-19_17-14-37/';
% Where is the behavior file
[Paths] = WhichComputer();
% where to look
[~,MouseName] = fileparts(fileparts(myKsDir));
[~,ephysfile] = fileparts(myKsDir);
Filename = strsplit(ephysfile,'_');
MatFileTag = regexprep(Filename{1},'-','');
if strcmp(MouseName(1:2),'NS')
    MouseName = MouseName(2:end);
end
rootfolder = fullfile(Paths.Grid.Behavior,MouseName);
BehaviorFiles = dir ([rootfolder,'/',[MouseName,'_',MatFileTag,'_r','*']]);
TuningFiles = dir ([rootfolder,'/',[MouseName,'_',MatFileTag,'_o','*']]);

% if ~isempty(BehaviorFiles) %& size(BehaviorFiles,1) == 1
%     % BehaviorPath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/Q9/Q9_20221119_r0.mat';
%     for s = 1:size(BehaviorFiles,1)
%         BehaviorPaths{s} = fullfile(rootfolder,BehaviorFiles(s).name);
%     end
% end
% if ~isempty(TuningFiles) %& size(TuningFiles,1) == 1
%     % TuningPath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/Q9/Q9_20221119_o0.mat';
%     for s = 1:size(TuningFiles,1)
%         TuningPaths{s} = fullfile(rootfolder,TuningFiles(s).name);
%     end
% end

if size(BehaviorFiles,1) == 2
    % first get the second session
    load(fullfile(rootfolder,BehaviorFiles(2).name),'session_data');
    TrialTrace = session_data.trace(:,7);
    TrialBinary = TrialTrace;
    TrialBinary(TrialBinary>0) = 1;
    behavior_idx = find(diff([0; TrialBinary; 0]));
    behavior_idx = reshape(behavior_idx,2,[])';
    behavior_TS = session_data.timestamps(behavior_idx);
    behavior_TS(:,3) = diff(behavior_TS')'; % trial durations
    BehaviorTrials_multi{2}.TimeStamps = behavior_TS;
end

% Then start with the first behavior file
BehaviorPath = fullfile(rootfolder,BehaviorFiles(1).name);

%% get Trial OFF Timestamps from the behavior file
load(BehaviorPath,'session_data');
TrialTrace = session_data.trace(:,7);
TrialBinary = TrialTrace;
TrialBinary(TrialBinary>0) = 1;
behavior_idx = find(diff([0; TrialBinary; 0]));
behavior_idx = reshape(behavior_idx,2,[])';
behavior_TS = session_data.timestamps(behavior_idx);
behavior_TS(:,3) = diff(behavior_TS')'; % trial durations
Odor_Ids = sqrt(TrialTrace(behavior_idx(1:5,1)));
TuningTrials = [];

%% get Trial Off Timestamps from the ephys side
if size(BehaviorFiles,1) == 2
    BehaviorTrials_multi{1}.TimeStamps = behavior_TS;
    [TTLs] = GetOepsAuxChannels_Multi(myKsDir, BehaviorTrials_multi, TuningTrials, [], 'KiloSorted', 1);
else
    [TTLs] = GetOepsAuxChannels(myKsDir, behavior_TS, TuningTrials, [], 'KiloSorted', 1);
end

%% calculate timestamp adjustment between closed loop behavior and ephys
ephys = TTLs.Trial(1:size(behavior_TS,1),2);
behavior = behavior_TS(:,2);
myfit = fit(behavior,ephys-behavior,'poly1');
TimestampAdjust.ClosedLoop(2) = myfit.p2; 
TimestampAdjust.ClosedLoop(1) = myfit.p1;

%% process sniff trace and create a digitized trace
[SniffTS, RespTrace] = ReadThermistorData(BehaviorPath);
% convert timestamps to ephys timestamps
SniffTS(:,1:3) = SniffTS(:,1:3) + ...
                 SniffTS(:,1:3)*TimestampAdjust.ClosedLoop(1) + ...
                 TimestampAdjust.ClosedLoop(2);
RespTrace(:,1) = RespTrace(:,1) + ...
                 RespTrace(:,1)*TimestampAdjust.ClosedLoop(1) + ...
                 TimestampAdjust.ClosedLoop(2);

% create a digitized sniff trace
TraceTS = (0:0.002:(TTLs.Trial(end,2)+1)); %@2ms resolution, ephys timebase
DigitalSniffs = TraceTS*0;
LocationSniffs = TraceTS*0;

for n = 1:size(SniffTS,1)
    idx(1) = find(TraceTS>=SniffTS(n,1),1,'first');
    idx(2) = find(TraceTS> SniffTS(n,2),1,'first') - 1;
    DigitalSniffs(idx(1):idx(2)) = 1;
    location = SniffTS(n,4);
    LocationSniffs(idx(1):idx(2)) = location;
end

% Traces.Timestamps{1} = TraceTS;
% Traces.SniffsDigitized{1} = DigitalSniffs;
% Traces.SniffsLocationed{1} = LocationSniffs;

%% create traces for Manifold and odor ON-OFF
Manifold = TraceTS*0;
valveTS = TTLs.AirManifold;
for n = 1:size(valveTS,1)
    idx(1) = find(TraceTS>=valveTS(n,1),1,'first');
    idx(2) = find(TraceTS> valveTS(n,2),1,'first') - 1;
    Manifold(idx(1):idx(2)) = 1;
end

Odor = TraceTS*0;
for x = 1:3
    switch x
        case 1
            valveTS = TTLs.Odor1;
        case 2
            valveTS = TTLs.Odor2;
        case 3
            valveTS = TTLs.Odor3;
    end
    for n = 1:size(valveTS,1)
        idx(1) = find(TraceTS>=valveTS(n,1),1,'first');
        idx(2) = find(TraceTS> valveTS(n,2),1,'first') - 1;
        Odor(idx(1):idx(2)) = x;
    end
end

Air = TraceTS*0;
valveTS = TTLs.Air;
if isnan(valveTS(1))
    valveTS(1) = 0;
end
for n = 1:size(valveTS,1)
    idx(1) = find(TraceTS>=valveTS(n,1),1,'first');
    idx(2) = find(TraceTS> valveTS(n,2),1,'first') - 1;
    Air(idx(1):idx(2)) = 1;
end

% find periods when neither odor nor air valves are open
allOFF = (~Odor.*~Air)';
Odor(find(allOFF)) = -1;

Traces.Manifold{1} = Manifold';
Traces.Odor{1} = Odor';
Traces.Trial{1} = Odor*0; % dummy trace


%% get Trial OFF Timestamps from the tuning file
if ~isempty(TuningFiles) %& size(TuningFiles,1) == 1
    if size(BehaviorFiles,1) == 1
        TuningPaths = fullfile(rootfolder,TuningFiles(1).name);
    elseif (BehaviorFiles(2).datenum - TuningFiles(1).datenum) > 0
        % make sure the chronology is correct
        TuningPaths = fullfile(rootfolder,TuningFiles(1).name);
    else
        keyboard;
    end
end

load(TuningPath,'session_data');
TrialTrace = session_data.trace(:,7);
TrialBinary = TrialTrace;
TrialBinary(TrialBinary>0) = 1;
tuning_idx = find(diff([0; TrialBinary; 0]));
tuning_idx = reshape(tuning_idx,2,[])';
if tuning_idx(end) > size(session_data.timestamps,1)
    tuning_idx(end) = size(session_data.timestamps,1);
end
tuning_TS = session_data.timestamps(tuning_idx);
tuning_TS(:,3) = diff(tuning_TS')'; % trial durations

%% find the right set of trials to match
behavior_last = size(behavior_TS,1);
offset = 0;
tuning_idx = [nan nan];
max_offset = size(TTLs.Trial,1) - size(tuning_TS,1) - behavior_last + 1;
while offset <= max_offset
    n1 = behavior_last + offset + 1;
    n2 = n1 + size(tuning_TS,1) - 1 - 2;
    if mean(1000*abs([tuning_TS(2:end-1,3)- TTLs.Trial(n1:n2,3)])) < 1
        tuning_idx = [n1 n2];
        break;
    end
    offset = offset + 1;
end

if isnan(tuning_idx(1))
    keyboard; % couldnot match trials between ephys and tuning
    n1 = behavior_last + 1;
    n2 = n1 + size(tuning_TS,1) - 1;
    [TTLs.Trial(n1:n2,3) tuning_TS(1:end,3)]
    uptowhich = 12;
    ephys = TTLs.Trial(n1+(1:uptowhich),2);
    behavior = tuning_TS(1:uptowhich,2);
    myfit = fit(behavior,ephys-behavior,'poly1');
    TimestampAdjust.Passive(2) = myfit.p2;
    TimestampAdjust.Passive(1) = myfit.p1;
else
    % calculate timestamp adjustment for the passive part
    ephys = TTLs.Trial(tuning_idx(1)-1:tuning_idx(2)+1,2);
    behavior = tuning_TS(1:end,2);
    myfit = fit(behavior,ephys-behavior,'poly1');
    TimestampAdjust.Passive(2) = myfit.p2;
    TimestampAdjust.Passive(1) = myfit.p1;
end

%% process sniff trace and create a digitized trace for the passive part
[SniffTS_passive, RespTrace_passive] = ReadThermistorData(TuningPath);
% convert timestamps to ephys timestamps
SniffTS_passive(:,1:3) = SniffTS_passive(:,1:3) + ...
                 SniffTS_passive(:,1:3)*TimestampAdjust.Passive(1) + ...
                 TimestampAdjust.Passive(2);
RespTrace_passive(:,1) = RespTrace_passive(:,1) + ...
                 RespTrace_passive(:,1)*TimestampAdjust.Passive(1) + ...
                 TimestampAdjust.Passive(2);

for n = 1:size(SniffTS_passive,1)
    if SniffTS_passive(n,2) < TraceTS(end)
        idx(1) = find(TraceTS>=SniffTS_passive(n,1),1,'first');
        idx(2) = find(TraceTS> SniffTS_passive(n,2),1,'first') - 1;
        DigitalSniffs(idx(1):idx(2)) = 1;
        location = SniffTS_passive(n,4);
        LocationSniffs(idx(1):idx(2)) = location;
    end
end

Traces.Timestamps{1} = TraceTS';
Traces.SniffsDigitized{1} = DigitalSniffs';
Traces.SniffsLocationed{1} = LocationSniffs';

%% Get Sniffs sorted for analyzing sniff tuning
[AllSniffs, ColumnInfo] = GetAllSniffs(Traces, 'quickmode', 1);
% mark sessionphase
f = find(AllSniffs(:,1)>RespTrace(end,1),1,'first');
AllSniffs(f:end,8) = AllSniffs(f:end,8) + 1;

%% if there's another round of behavior stuff?
if size(BehaviorFiles,1) == 2
    keyboard;
end
%%
%end
