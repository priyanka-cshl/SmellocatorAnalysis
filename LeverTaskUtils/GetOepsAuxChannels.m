function [TTLs,ReplayTTLs,EphysTuningTrials,AuxData] = GetOepsAuxChannels(myKsDir, BehaviorTrials, TuningTrials, TrialSequence, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('ADC', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('KiloSorted', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
GetAux = params.Results.ADC;
UseSortingFolder = params.Results.KiloSorted;     

%% defaults
global SampleRate; % Behavior Acquisition rate

%% Get Trial Timestamps from the OpenEphys Events file
if UseSortingFolder
    if exist(fullfile(myKsDir,'myTTLfile_1.mat'))
        foo             = load(fullfile(myKsDir,'myTTLfile_1.mat'),'TTLs');
        data            = foo.TTLs.data;
        timestamps      = foo.TTLs.timestamps;
        info.eventId    = foo.TTLs.info.eventId;
        offset          = foo.TTLs.offset;
    else
        TTLs = [];
        EphysTuningTrials = [];
        AuxData = [];
        ReplayTTLs = [];
        disp('empty events file');
        return
    end
else
    filename = fullfile(myKsDir,'all_channels.events');

    % hack to avoid going through empty files
    temp = dir(filename);
    if ~isempty(temp) % old GUI
        if ~temp.bytes
            TTLs = [];
            EphysTuningTrials = [];
            AuxData = [];
            ReplayTTLs = [];
            disp('empty events file');
            return
        end

        [data, timestamps, info] = load_open_ephys_data(filename); % data has channel IDs

        % adjust for clock offset between open ephys and kilosort
        [offset] = AdjustClockOffset(myKsDir);
    else % new GUI
        session = Session(fileparts(myKsDir));
        Events = session.recordNodes{1}.recordings{1}.ttlEvents('Acquisition_Board-100.Rhythm Data');
        data            = Events.channel;
        timestamps      = Events.timestamp;
        info.eventId    = Events.state;

        % adjust for clock offset between open ephys and kilosort
        offset = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps(1);
    end
end

timestamps = timestamps - offset;

%% Get various events
TTLTypes = unique(data);
Tags = {'Air', 'Odor1', 'Odor2', 'Odor3', 'Trial', 'Reward', 'AirManifold', 'Licks'};
for i = 1:numel(TTLTypes)
    On = timestamps(intersect(find(info.eventId),find(data==TTLTypes(i))));
    Off = timestamps(intersect(find(~info.eventId),find(data==TTLTypes(i))));
    % delete the first off value, if it preceeds the On
    Off(Off<On(1)) = [];
    On(On>Off(end)) = [];
    
    if length(On)>length(Off)
        keyboard;
        foo = [On(1:end-1) Off]';
        goo = foo(:);
        On(find(diff(goo)<0,1,'first')/2,:) = [];
    end
        
    temp = [On Off Off-On];
    
    % ignore any transitions faster than 1 ms - behavior resolution is 2 ms
    temp(temp(:,3)<0.001,:) = [];
    
    % sometimes lines can toggle high low quickly leading to
    % splitting of one trial into two
    splitTrials = find(abs(Off(:,1) - circshift(On(:,1),-1))<0.001);
    if any(splitTrials)
        disp(['merging ',num2str(numel(splitTrials)),' split ',char(Tags(i)),' Trials']);
        Off(splitTrials,:)  = [];
        On(splitTrials+1,:) = [];
    end
    
    temp = [On Off Off-On];
    TTLs.(char(Tags(i))) = temp;
end

if size(TTLs.Trial,1)<5
    TTLs = [];
    EphysTuningTrials = [];
    AuxData = [];
    disp('no valid oeps file found');
    return
end

%% mismatch between behavior and oeps trials on first trial
if ~isempty(BehaviorTrials)
    if any(abs(BehaviorTrials(1:5,3)-TTLs.Trial(1:5,3))>=0.01)
        % case 1 : behavior file has an extra trial
        if ~any(abs(BehaviorTrials(2:6,3)-TTLs.Trial(1:5,3))>=0.01)
            % append an extra NaN trial on the Oeps side
            TTLs.Trial = [NaN*TTLs.Trial(1,:); TTLs.Trial];
            
            % case 2 : behavior acq started mid-trial, first trial might be a bit shorter
        elseif ~any(abs(BehaviorTrials(2:6,3)-TTLs.Trial(2:6,3))>=0.01)
            % do nothing
        end
    end
    
    % Is it the right recording file?
    if size(TTLs.Trial,1)<size(BehaviorTrials,1)
        TTLs = [];
        EphysTuningTrials = [];
        AuxData = [];
        disp('behavior and ephys files do not match');
        return
    end
    
    y = corrcoef(BehaviorTrials(2:20,3),TTLs.Trial(2:20,3));
    if y(2,1)<0.99
        TTLs = [];
        EphysTuningTrials = [];
        AuxData = [];
        disp('behavior and ephys files do not match');
        return
    end
end

%% find the odor ON time 
count = 0;
ReplayTTLs = [];
for i = 1:size(TTLs.Trial,1) % every trial
    
     % find the last odor valve ON transition just before this trial start
     if i > 1
         t1 = TTLs.Trial(i-1,2);
     else
         t1 = 0;
     end
     t2 = TTLs.Trial(i,1);
     
     ValveEvents = [];
     for thisOdor = 1:3
         myEvents = intersect(find(TTLs.(['Odor',num2str(thisOdor)])(:,1)>t1),...
             find(TTLs.(['Odor',num2str(thisOdor)])(:,1)<t2));
         myTimeStamps = TTLs.(['Odor',num2str(thisOdor)])(myEvents,:);
         ValveEvents = vertcat(ValveEvents,...
             [myTimeStamps thisOdor*ones(numel(myEvents),1)]);
     end
     
     if ~isempty(ValveEvents)         
         [t3,x] = max(ValveEvents(:,1));
         TTLs.Trial(i,4:5) = [t3-t2 ValveEvents(x,4)];
     else
         TTLs.Trial(i,4:5) = [NaN 0];
     end
             
    % for replay trials
    %if TTLs.Trial(i,3) > 1 + mode(BehaviorTrials(:,3))
    if ~isempty(BehaviorTrials)
        if TTLs.Trial(i,3) > 5*mean(BehaviorTrials(:,3)) % at least 5 trials in the replay
            tstart = TTLs.Trial(i,1);
            tstop  = TTLs.Trial(i,2);
            O1 = intersect(find(TTLs.Odor1(:,1)>tstart),find(TTLs.Odor1(:,1)<tstop));
            O2 = intersect(find(TTLs.Odor2(:,1)>tstart),find(TTLs.Odor2(:,1)<tstop));
            O3 = intersect(find(TTLs.Odor3(:,1)>tstart),find(TTLs.Odor3(:,1)<tstop));
            
            % keep the odor transition that happened just before trial start
            if ~isempty(ValveEvents)
                ValveEvents = ValveEvents(x,:);
            end
            for j = 1:3
                temp = eval(['TTLs.Odor',num2str(j),'(O',num2str(j),',:);']);
                ValveEvents = vertcat(ValveEvents,...
                    [temp j*ones(size(temp,1),1)]);
            end
            % reference odor transitions w.r.t. trial ON
            ValveEvents(:,1:2) = ValveEvents(:,1:2) - t2 ;
            [~,sortID] = sort(ValveEvents(:,1));
            
            if size(ValveEvents,1)>5
                count = count + 1;
                ReplayTTLs.OdorValve{count} = ValveEvents(sortID,:);
                ReplayTTLs.TrialID(count)   = i;
            end
            
        end
    end
    
end

% Align Passive Tuning trials
if ~isempty(TuningTrials)
    EphysTuningTrials = AlignPassiveTuningTrials(TuningTrials, TTLs, size(BehaviorTrials,1), TrialSequence);
else
    EphysTuningTrials = [];
end

AuxData = [];
if GetAux
    if UseSortingFolder
        %if exist(fullfile(myKsDir,'myTTLfile_1.mat'))
        
    else

        %% Get analog/digital AuxData from Oeps files - for comparison with behavior data
        foo = dir(fullfile(myKsDir,'*_ADC1.continuous')); % pressure sensor
        filename = fullfile(myKsDir,foo.name);
        [Auxdata1, timestamps, ~] = load_open_ephys_data(filename); % data has channel IDs
        foo = dir(fullfile(myKsDir,'*_ADC2.continuous')); % thermistor
        filename = fullfile(myKsDir,foo.name);
        [Auxdata2, ~, ~] = load_open_ephys_data(filename); % data has channel IDs

        % adjust for clock offset between open ephys and kilosort
        timestamps = timestamps - offset;

        % downsample to behavior resolution

        AuxData(:,1) = 0:1/SampleRate:max(timestamps);
        AuxData(:,2) = interp1q(timestamps,Auxdata1,AuxData(:,1)); % pressure sensor
        AuxData(:,3) = interp1q(timestamps,Auxdata2,AuxData(:,1)); % thermistor
        % create a continuous TrialOn vector
        for MyTrial = 1:size(TTLs.Trial,1)
            [~,start_idx] = min(abs(AuxData(:,1)-TTLs.Trial(MyTrial,1)));
            [~,stop_idx]  = min(abs(AuxData(:,1)-TTLs.Trial(MyTrial,2)));
            AuxData(start_idx:stop_idx,4) = 1;
        end
    end
end
end
