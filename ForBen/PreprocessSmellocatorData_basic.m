%% function to parse behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

function [] = PreprocessSmellocatorData_basic(MyFilePath,overwriteflag)

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
Paths.Grid.Behavior_processed   = '/home/priyanka/Dropbox/Smellocator_Behavior'; %'/mnt/grid-hs/pgupta/Smellocator/Behavior';

savepath = fullfile(Paths.Grid.Behavior_processed,AnimalName,[MyFileName,'_processed.mat']);
if ~overwriteflag && exist(savepath)
    reply = input('This session has already been processed. \nDo you want to overwrite? Y/N [Y]: ','s');
    if strcmp(reply,'N')
        return;
    end
end

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
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
Traces.OdorLocation = Traces.Motor;
Traces.TrialState = Traces.Trial;
Traces = rmfield(Traces,{'Motor'; 'Trial'});

%     {'Lever'       }
%     {'Sniffs'      }
%     {'Licks'       }
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

save(savepath, 'Traces', 'TrialInfo', 'TargetZones', 'startoffset', 'SampleRate');
    
end