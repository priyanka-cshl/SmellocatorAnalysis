
function [TracesOut,SingleUnits] = LoadSniffSessions_wdw(WhereSession,savemode)

if nargin < 2
    savemode = 0;
end

%WhereSession = '/mnt/data/Sorted/T3/2025-05-08_16-20-35/quickprocesssniffs.mat';

%% Load the relevant variables
% Sniffs
[SniffCoords,SniffProps] = TallyThermistorNPressureSniffs(WhereSession);
if ~isempty(SniffCoords)
    load(WhereSession, 'RespirationData');
else
    load(WhereSession, 'CuratedMFSSniffTimestamps','RespirationData');
    CuratedSniffTimestamps = CuratedMFSSniffTimestamps;
end

try
    load(WhereSession,'quickprocessOdorTTLs.mat');
    if ~exist('TTLs','var')
        TTLs = [];
    end
catch
    TTLs = [];
end

myKsDir = fileparts(WhereSession);
SingleUnits = GetSingleUnits(myKsDir, 3);

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
end
% fake trace for fitting kernels similar to that for the lever task
TracesOut.Motor{1} = TracesOut.Timestamps{1}*0;

% sniffing specific
% add a filtered sniff trace
TracesOut.SniffsFiltered{1}     = RespirationData(:,3);
TracesOut.MassFlowSensor{1}     = RespirationData(:,4) - 2.5;
TracesOut.SniffsLocationed{1}   = TracesOut.SniffsFiltered{1}*nan;

if ~isempty(SniffCoords)
    % make several digitized sniffs
    DigitalSniffs = zeros(numel(TracesOut.SniffsFiltered{1}),3);
    for n = 1:size(SniffCoords,1)
        idx = SniffCoords(n,[4 5 9 10 14 15]);
        idx = reshape(idx,2,[])';
        for m = 1:3
            DigitalSniffs(idx(m,1):idx(m,2),m) = 1;
        end
    end
    TracesOut.SniffsDigitized{1} = DigitalSniffs(:,1);
    TracesOut.MFS2ThermDigitized{1} = DigitalSniffs(:,2);
    TracesOut.MFSDigitized{1} = DigitalSniffs(:,3);
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

%%
if savemode
    [~,MouseName] = fileparts(fileparts(myKsDir));
    [~,filename] = fileparts(myKsDir);
    filename = [MouseName,'_',regexprep(filename(1,1:10),'-',''),'_r0_processed.mat']; 
    savepath = '/home/priyanka/Desktop';
    save(fullfile(savepath,'forWDW',filename),'TracesOut','SingleUnits');
end

end
