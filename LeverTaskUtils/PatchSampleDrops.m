%function [] = PatchSampleDrops(MyFilePath)
%% function to patch samples in the analog traces that got dropped
% from the parallel OEPS recordings

Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-matlab-tools'))); % for the new OEPS GUI: https://github.com/open-ephys/open-ephys-matlab-tools (use commit 10-04-22)

if exist(MyFilePath) & ~isempty(fileparts(MyFilePath))
    [FilePaths, MyFileName] = fileparts(MyFilePath);
    [~,AnimalName] = fileparts(FilePaths);
else
    foo = regexp(MyFilePath,'_','split');
    AnimalName = foo{1};
    MyFilePath = fullfile(Paths.Grid.Behavior,AnimalName,MyFilePath);
    [FilePaths, MyFileName] = fileparts(MyFilePath); %#ok<*ASGLU>
end

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
drop_points = find(abs(diff(MyData(:,1)))>0.01); % indices at which timestamps were dropped

if ~any(MyData(drop_points+1,6)==0) && ~any(MyData(drop_points,6)~=0)
    % all timestamp drops are st trial starts
    
    % find the corresponding ephys file
    [myEphysdir] = WhereRawOEPS(MyFileName,FilePaths);
    
    %% read the open ephys flat binary file
    session = Session(myEphysdir);
    TotalSamples = size(session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps,1);
    [TTLs] = OEPS_getTTLs_newGUI(session);
    
    TS = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps;
    TS = TS - session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps(1);
        
    % first find the timestamp adjustment that needs to be made
    TrialOn.behavior = MyData(find(diff(MyData(:,6))>0,5,'first'),1);
    
    if ~any((diff(TrialOn.behavior) - diff(TTLs.Trial(1:5,1))) > 0.05)
        TimestampAdjust.ClosedLoop = median(TTLs.Trial(1:5,1) - TrialOn.behavior);
    else
        keyboard;
    end
    
    % Lever, OdorLocation, Sniffs, TrialState, Rewards, Licks, Timestamps
    LeverOEPS = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(43,:);
    SniffsOEPS
    % rescale 
    LeverOEPS = 2.5 + (double(LeverOEPS)/12800);
    
    %%
    drop_points = find(abs(diff(MyData(:,1)))>0.01, 1, 'first'); % first indices at which timestamps were dropped
    while ~isempty(drop_points)
        idx_B = [drop_points drop_points+1 (find(MyData((drop_points+1):end,6)==0,1,'first') + drop_points)];
        ts_B  = MyData(idx_B,1);
        ts_O  = ts_B + TimestampAdjust.ClosedLoop;
        % find the trace indices in OEPS that are closest to this
        for i = 1:3
            [~, idx_O(i)] = min(abs(TS-ts_O(i)));
        end
        
        % interpolate the OEPS trace 
        Xq = MyData(idx_B(2):idx_B(3),1);
        X  = TS(idx_O(2):idx_O(3));
        V  = LeverOEPS(idx_O(2):idx_O(3));
        Vq = interp1(X-X(1),V,Xq-Xq(1));
        
        LeverSnippet_B  = MyData(idx_B(2):idx_B(3),4);
        % account for trace_offset
        trace_offset = nanmedian(Vq-LeverSnippet_B);
        LeverSnippet_O  = Vq - trace_offset;
        
        if abs(nanmedian(LeverSnippet_O - LeverSnippet_B)) < 0.02
            % interpolate again
            X2q = MyData(idx_B(1),1):0.002:MyData(idx_B(2),1);
            X2  = TS(idx_O(1):idx_O(2));
            V2  = LeverOEPS(idx_O(1):idx_O(2));
            V2q = interp1(X2-X2(1),V2,X2q-X2q(1)) - trace_offset;
            
            MissingTS = X2q(2:end-1);
            MissingSnippet = V2q(2:end-1);
            MissingBlock = nan(numel(MissingTS),size(MyData,2));
            MissingBlock(:,1) = MissingTS;
            MissingBlock(:,4) = MissingSnippet;
            MissingBlock(:,6) = MyData(drop_points+1,6);
            
            MyData = vertcat(MyData(1:drop_points,:), ...
                            MissingBlock, ...
                            MyData((drop_points+1):end,:) ); 
        end
        
        drop_points = find(abs(diff(MyData(:,1)))>0.01, 1, 'first'); % first indices at which timestamps were dropped
    end
    
    %LeverSnippet 
end

function [MyData] = patch_smellocatordata_gap(drop_point,MyData)
% drop_points = last index in MyData after which the samples got dropped.
% MyData = actual old data with the timestamp drops

% 1. create a dummy block
ts_before   = MyData(drop_point,1);
ts_after    = MyData(drop_point+1,1);
ts_dropped  = (ts_before + 0.002):0.002:(ts_after - 0.002);
DummyBlock  = nan(numel(ts_dropped),size(MyData,2));

% 2. find trial off indices and timestamps
idx_trialoff    = find(MyData((drop_point+1):end,6)==0,1,'first') + drop_point;
ts_full         = [ts_dropped MyData((drop_point+1):idx_trialoff,1)'];

% 3. convert to OEPS indices to get the correct trace portions out
[~,idx_OEPS(1)] = min(abs(TS-ts_full(1)-TimestampAdjust.ClosedLoop));
[~,idx_OEPS(2)] = min(abs(TS-ts_full(end)-TimestampAdjust.ClosedLoop));
idx_ephys = idx_OEPS(1):idx_OEPS(2);
ts_OEPS   = TS(idx_ephys);
% Thermistor, Lever, Piezo, RE
Traces_OEPS = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').samples(42:45,:);

% 4. interpolate to get traces at behavior resolution
% Vq = interp1(X,V,Xq)
Traces_interp = interp1( (ts_OEPS-ts_OEPS(1)), (Traces_OEPS), (ts_full-ts_full(1)) );

% 5. get behavior traces to confirm the match
whichcolumns = [15 4 10 5]; % Thermistor, Lever, Piezo, RE
Data_B = MyData((drop_point+1):idx_trialoff,whichcolumns);
% which portion of the interpolated traces to compare?
[~,idx_split] = min(abs(ts_full-ts_after));
Data_E = Traces_interp(idx_split:end,:);

% 6. compare

end
%end