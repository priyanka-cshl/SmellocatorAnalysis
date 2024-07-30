%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

% same as PreprocessSmellocatorData_basic, except with sampledrops patched
% in from the OEPS file when possible

function [] = PreprocessSmellocatorData_basic_v2(MyFilePath,overwriteflag)

if nargin<2
    overwriteflag = 0;
end

%% Add relevant repositories
Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));
addpath(genpath(fullfile(Paths.Code,'open-ephys-matlab-tools'))); % for the new OEPS GUI: https://github.com/open-ephys/open-ephys-matlab-tools (use commit 10-04-22)
addpath(genpath(fullfile(Paths.Code,'MatlabUtils')));
addpath(genpath(fullfile(Paths.Code,'npy-matlab/'))); % path to npy-matlab scripts

%% globals
% global MyFileName;
global SampleRate;
SampleRate = 500; % Samples/second
global startoffset;
startoffset = 1; % in seconds
global savereplayfigs;
savereplayfigs = 0;
global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips]
errorflags = [0 0 0 0];
global TargetZones; %#ok<*NUSED>

%% core data extraction (and settings)
if ~exist(MyFilePath)
    foo = regexp(MyFilePath,'_','split');
    AnimalName = foo{1};
    MyFilePath = fullfile(Paths.Grid.Behavior,AnimalName,MyFilePath);
    [FilePaths, MyFileName] = fileparts(MyFilePath); %#ok<*ASGLU>
else
    [FilePaths, MyFileName] = fileparts(MyFilePath);
    [~,AnimalName] = fileparts(FilePaths);
end

%% check if the preprocessed version already exists - locally or on the server
Paths.Grid.Behavior_processed   = '/home/priyanka/Dropbox/Smellocator_Behavior_2';
%Paths.Grid.Behavior_processed   = '/mnt/data/Processed/Behavior'; %'/home/priyanka/Dropbox/Smellocator_Behavior'; %'/mnt/grid-hs/pgupta/Smellocator/Behavior';

savepath = fullfile(Paths.Grid.Behavior_processed,AnimalName,[MyFileName,'_processed.mat']);
if ~overwriteflag && exist(savepath)
    reply = input('This session has already been processed. \nDo you want to overwrite? Y/N [Y]: ','s');
    if strcmp(reply,'N')
        return;
    end
end

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath); %,'PIDmode',2); % 2 for populating lick column with lick piezo instead of binary licks
FileLocations.Behavior = MyFilePath;
[FilePaths, MyFileName] = fileparts(MyFilePath);
disp(MyFileName);

%% Parse into trials
%[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[Trials,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

% sanity check - did some guess work in CorrectMatlabSampleDrops to compute
% odor start - check if it made sense
if ~isempty(InitiationsFixed)
    if any(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01)
        weirdo = find(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01);
        if any(TrialInfo.OdorStart(InitiationsFixed(weirdo),2)>-1)
            disp('something funky with computing odorstart from Lever trace');
            %keyboard;
            TrialInfo.OdorStart(InitiationsFixed(weirdo),2) = TrialInfo.OdorStart(InitiationsFixed(weirdo),1);
        else
            % Initiation hold was larger than a second - that's couldn't
            % compute it accurately from trial traces
            TrialInfo.OdorStart(InitiationsFixed(weirdo),2) = TrialInfo.OdorStart(InitiationsFixed(weirdo),1);
        end
    end
end

%% Saving stuff in one place
if ~exist(fileparts(savepath),'dir')
    mkdir(fileparts(savepath));
end

%% For Ben and Ryan - remove unnecessary stuff
TrialInfo.OdorStart(:,2) = [];
TrialInfo.Success(:,2) = [];
Traces.OdorLocation     = Traces.Motor;
Traces.TrialState       = Traces.Trial;
Traces.LicksPiezo       = Traces.Piezo;
Traces.LicksBinary      = Traces.Licks;
Traces = rmfield(Traces,{'Motor'; 'Trial'; 'Licks'; 'Piezo'});

%     {'Lever'       }
%     {'Sniffs'      }
%     {'LicksPiezo'  }
%     {'LicksBinary' }
%     {'Rewards'     }
%     {'Timestamps'  }
%     {'OdorLocation'}
%     {'TrialState'  }

extrafields = {'Offset'; 'TraceIndices'; 'TraceDuration'; 'SessionIndices'; 'SessionTimestamps'; ...
    'TimeIndices'; 'Timestamps'; 'Duration'; 'Reward'; 'TransferFunctionLeft'};
TrialInfo = rmfield(TrialInfo,extrafields);
TrialInfo.TrialID = TrialInfo.TrialID';

%     {'TrialID'             }
%     {'TimeStampsDropped'   }
%     {'Odor'                }
%     {'OdorStart'           }
%     {'TargetZoneType'      }
%     {'Success'             }
%     {'InZone'              }
%     {'HoldSettings'        }
%     {'Perturbation'        }

%% sniff peak-valley detection
SniffTS = [];
for trialID = 1:numel(TrialInfo.Odor)
    % detect sniff timestamps per trial
    TraceTimeStamps = Traces.Timestamps{trialID};
    Thermistor      = Traces.Sniffs{trialID};
    OdorLocation    = Traces.OdorLocation{trialID};
    thisTrialSniffs = ProcessThermistorData([TraceTimeStamps Thermistor OdorLocation]);

    % apppend to AllSniffs and also handle end cases
    if trialID > 1
        % if first sniff started earlier than this trace
        if isnan(thisTrialSniffs(1,1)) 
            [~,f] = min(sum(abs(SniffTS(:,2:3) - thisTrialSniffs(1,2:3)),2));
            if ~isempty(f)
                thisTrialSniffs(1,:) = SniffTS(f,1:4);
            end
        end
        % if the last sniff in Allsniffs was truncated
        if isnan(SniffTS(end,3))
            [~,f] = min(sum(abs(thisTrialSniffs(:,1:2) - SniffTS(end,1:2)),2));
            if ~isempty(f)
                SniffTS(end,1:4) = thisTrialSniffs(f,1:4);
            end
        end
        % flag overlapping sniffs
        f = find(ismember(SniffTS(:,1:4),thisTrialSniffs(1,1:4),'rows'));
        if ~isempty(f)
            SniffTS(f:end,5) = -SniffTS(f:end,5);
        end
    end

    SniffTS = vertcat(SniffTS, ...
                        [thisTrialSniffs (trialID + 0*thisTrialSniffs(:,1))]);
end

save(savepath, 'Traces', 'TrialInfo', 'TargetZones', 'SniffTS', 'startoffset', 'SampleRate');
    
end