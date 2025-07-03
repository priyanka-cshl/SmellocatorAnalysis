function [AllSniffs] = QuickSniffTTLMapper_v2(myKsDir)
% myKsDir = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Ephys/Q9/2022-11-19_17-14-37/';
% location where the sorted output is - contains cluster info and a TTL
% file generated during Kilosort

%% 1: Where is the corresponding behavior/tuning file
%   - need this to load trials to calculate timestamp difference
%   - and then process respiration with corrected timestamps

% where to look
[Paths] = WhichComputer();
if strcmp(myKsDir(end),filesep)
    myKsDir = myKsDir(1:end-1);
end

% have sniffs already been processed
if exist(fullfile(myKsDir,'quickprocesssniffs.mat'))
    load (fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs');
    if strcmp(myKsDir,'/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31')
        AllSniffs(find(AllSniffs(:,1)>9000),:) = [];
    end
    return;
end

[~,ephysfile] = fileparts(myKsDir);
[~,MouseName] = fileparts(fileparts(myKsDir));
Filename = strsplit(ephysfile,'_'); % splits date and timestamp
MatFileTag = regexprep(Filename{1},'-',''); % removes dashes from date
% hack for dealing with resorting of sessions from batch S from Marie
if strcmp(MouseName(1:2),'NS')
    MouseName = MouseName(2:end);
end
rootfolder = fullfile(Paths.Grid.Behavior,MouseName); % where behavior and tuning files are saved
BehaviorFiles = dir ([rootfolder,'/',[MouseName,'_',MatFileTag,'_r','*']]);
TuningFiles = dir ([rootfolder,'/',[MouseName,'_',MatFileTag,'_o','*']]);


%% for every behavior (Closed loop) session
% get Trial OFF Timestamps from the behavior file
if size(BehaviorFiles,1) == 2
    % first get the second session
    BehaviorPath = fullfile(rootfolder,BehaviorFiles(2).name);
    BehaviorTrials_multi{2}.TimeStamps = TrialsfromMatFile(BehaviorPath);

    % also check chronology
    if ~isempty(TuningFiles)
        if (BehaviorFiles(1).datenum - TuningFiles(1).datenum < 0) && ...
            (BehaviorFiles(2).datenum - TuningFiles(1).datenum > 0)
            if size(TuningFiles,1) > 1
                if (BehaviorFiles(2).datenum - TuningFiles(2).datenum > 0)
                    disp('File order does not make sense');
                    keyboard;
                end
            end
        else
            disp('File order does not make sense');
            keyboard;
        end
        
         
    end
end

% start with the first session
BehaviorPath = fullfile(rootfolder,BehaviorFiles(1).name);
[behavior_TS] = TrialsfromMatFile(BehaviorPath);

%% get Trial Off Timestamps from the ephys side
if size(BehaviorFiles,1) == 2
    BehaviorTrials_multi{1}.TimeStamps = behavior_TS;
    [TTLs] = GetOepsAuxChannels_Multi(myKsDir, BehaviorTrials_multi, [], [], 'KiloSorted', 1);
else
    [TTLs] = GetOepsAuxChannels(myKsDir, behavior_TS, [], [], 'KiloSorted', 1);
end

%% create a timestamps trace @2ms resolution, ephys timebase to use for
% creating other behavior related traces
TraceTS = (0:0.002:(TTLs.Trial(end,2)+1));
% also create a blank  digitized and locationed sniff trace
DigitalSniffs = TraceTS*0;
LocationSniffs = TraceTS*0;

%% create traces for Manifold and odor ON-OFF
Manifold = TraceTS*0;
ValveTS = TTLs.AirManifold;
[Manifold] = TStoValveTrace(ValveTS,Manifold,1,TraceTS);

Odor = TraceTS*0;
for x = 1:3
    ValveTS = TTLs.(['Odor',num2str(x)]);
    [Odor] = TStoValveTrace(ValveTS,Odor,x,TraceTS);
end

Air = TraceTS*0;
ValveTS = TTLs.Air;
if isnan(ValveTS(1)) % hack
    ValveTS(1) = 0;
end
[Air] = TStoValveTrace(ValveTS,Air,1,TraceTS);

% find periods when neither odor nor air valves are open
allOFF = (~Odor.*~Air)';
Odor(find(allOFF)) = -1;

Traces.Manifold{1} = Manifold';
Traces.Odor{1} = Odor';
Traces.Trial{1} = Odor*0; % dummy trace (not used)

%% calculate timestamp adjustment between closed loop behavior and ephys
start_from = 1; % first session so start looking from the beginning of the trial list
[TSadjust] = matchtrialsandclockdiff(behavior_TS, TTLs.Trial, start_from);

%% process the respiration data from the matlab file, detect peaks/valleys, convert TS to OEPS and return digitized and locationed sniff traces
[DigitalSniffs, LocationSniffs, SessMarker(1)] = SniffProcess(BehaviorPath,TSadjust,TraceTS,DigitalSniffs,LocationSniffs);

%% get Trials from the matlab tuning file to calculate the TS correction for tuning sniffs
if ~isempty(TuningFiles) %& size(TuningFiles,1) == 1
    TuningPath = fullfile(rootfolder,TuningFiles(1).name);
    [tuning_TS] = TrialsfromMatFile(TuningPath);

    %% find the right set of trials to match
    start_from = size(behavior_TS,1);
    [TSadjust, trialsdone] = matchtrialsandclockdiff(tuning_TS, TTLs.Trial, start_from);

    %% process the respiration data from the matlab file, detect peaks/valleys, convert TS to OEPS and return digitized and locationed sniff traces
    [DigitalSniffs, LocationSniffs, SessMarker(2)] = SniffProcess(TuningPath,TSadjust,TraceTS,DigitalSniffs,LocationSniffs);
end

%% if there's another round of behavior stuff?
if size(BehaviorFiles,1) == 2
    behavior_TS = BehaviorTrials_multi{2}.TimeStamps;
    [TSadjust, trialsdone] = matchtrialsandclockdiff(behavior_TS, TTLs.Trial, trialsdone);

    % process respiration
    BehaviorPath = fullfile(rootfolder,BehaviorFiles(2).name);
    [DigitalSniffs, LocationSniffs, SessMarker(3)] = SniffProcess(BehaviorPath,TSadjust,TraceTS,DigitalSniffs,LocationSniffs);

    %% another round of tuning as well?
    if size(TuningFiles,1) == 2
        TuningPath = fullfile(rootfolder,TuningFiles(2).name);
        [tuning_TS] = TrialsfromMatFile(TuningPath);
        [TSadjust, trialsdone] = matchtrialsandclockdiff(tuning_TS, TTLs.Trial, trialsdone);
        [DigitalSniffs, LocationSniffs, SessMarker(4)] = SniffProcess(TuningPath,TSadjust,TraceTS,DigitalSniffs,LocationSniffs);

    end
end

%%
Traces.Timestamps{1} = TraceTS';
Traces.SniffsDigitized{1} = DigitalSniffs';
Traces.SniffsLocationed{1} = LocationSniffs';

%% Get Sniffs sorted for analyzing sniff tuning
[AllSniffs, ColumnInfo] = GetAllSniffs(Traces, 'quickmode', 1);
% mark sessionphase
for s = 1:numel(SessMarker)-1
    f = find(AllSniffs(:,1)>SessMarker(s),1,'first');
    AllSniffs(f:end,8) = AllSniffs(f:end,8) + 1;
end

save(fullfile(myKsDir,'quickprocesssniffs.mat'),'AllSniffs');

%% function definitions

    function [trial_TS] = TrialsfromMatFile(Filepath)
        % get Trial OFF Timestamps from the behavior/tuning file
        load(Filepath,'session_data');

        % load and binzarize the trial ON-OFF trace
        MyTrialTrace = session_data.trace(:,7);
        MyTrialBinary = MyTrialTrace;
        if any(MyTrialBinary<0)
            MyTrialBinary(MyTrialBinary~=0) = 1;
        else
            MyTrialBinary(MyTrialBinary>0) = 1;
        end
        % detect trial transitions and get transition timestamps in behavior space
        trial_idx = find(diff([0; MyTrialBinary; 0]));
        trial_idx = reshape(trial_idx,2,[])';
        if trial_idx(end) > size(session_data.timestamps,1)
            trial_idx(end) = size(session_data.timestamps,1);
        end
        trial_TS = session_data.timestamps(trial_idx);
        trial_TS(:,3) = diff(trial_TS')'; % trial durations
    end

    function [TSadjust, trialsdone] = matchtrialsandclockdiff(Trials_B, Trials_E, start_E)
        % ephys trials can be superset - so find the correct matching set
        offset = 0;
        match_idx = [nan nan];
        max_offset = 5;
    
        while offset <= max_offset
            n1 = start_E + offset; % starting index in the ephys list
            n2 = n1 + size(Trials_B,1) - 1; % last index in the ephys list
            if n2 > size(Trials_E,1)
                break;
            end
            if mean(1000*abs([Trials_B(2:end-1,3)- Trials_E((n1+1):(n2-1),3)])) < 1 % less than 1 ms difference in durations
                match_idx = [n1 n2];
                break;
            end
            offset = offset + 1;
        end
    
        if ~isnan(match_idx(1))
            ephys = Trials_E(match_idx(1):match_idx(2),2);
            behavior = Trials_B(:,2);
            trialsdone = match_idx(2);
        else
            disp('could not match OEPS and matlab trial lists');
            keyboard;
            n1 = start_E + 1;
            n2 = n1 + size(Trials_B,1) - 1;
            [Trials_E(n1:n2,3) Trials_B(1:end,3)]
            uptowhich = 12;
            [Trials_E(n1+(1:uptowhich),3) Trials_B((1:uptowhich),3)]
            ephys = Trials_E(n1+(1:uptowhich),2);
            behavior = Trials_B(1:uptowhich,2);
            trialsdone = n2;
        end
        myfit = fit(behavior,ephys-behavior,'poly1');
        TSadjust(2) = myfit.p2;
        TSadjust(1) = myfit.p1;
    end
    
    function [ValveTrace] = TStoValveTrace(ValveTS,ValveTrace,ONvalue,TraceTS) 
        for v = 1:size(ValveTS,1)
            vx(1) = find(TraceTS>=ValveTS(v,1),1,'first');
            vx(2) = find(TraceTS> ValveTS(v,2),1,'first') - 1;
            ValveTrace(vx(1):vx(2)) = ONvalue;
        end
    end

    function [DigitalSniffs, LocationSniffs, SessionLength] = SniffProcess(BehaviorPath,TSadjust,TraceTS,DigitalSniffs,LocationSniffs)
        % process sniff trace and create a digitized trace
        [SniffTS, RespTrace] = ReadThermistorData(BehaviorPath);
        % convert timestamps to ephys timestamps
        SniffTS(:,1:3) = SniffTS(:,1:3) + SniffTS(:,1:3)*TSadjust(1) + TSadjust(2); % convert timestamps to OEPS
        RespTrace(:,1) = RespTrace(:,1) + RespTrace(:,1)*TSadjust(1) + TSadjust(2); % convert timestamps to OEPS
        SessionLength = RespTrace(end,1);

        for n = 1:size(SniffTS,1)
            idx(1) = find(TraceTS>=SniffTS(n,1),1,'first');
            idx(2) = find(TraceTS> SniffTS(n,2),1,'first') - 1;
            DigitalSniffs(idx(1):idx(2)) = 1;
            location = SniffTS(n,4);
            LocationSniffs(idx(1):idx(2)) = location;
        end
    end


end
