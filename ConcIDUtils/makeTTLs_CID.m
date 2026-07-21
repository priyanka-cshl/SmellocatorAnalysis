function [TTLs,StimSettings] = makeTTLs_CID(myKsDir,myStimFile) %(myDir,myStimFile)

% myKsDir = '/mnt/data/Sorted/T2/_2025-05-21_09-18-56_2025-05-21_11-11-12_2025-05-21_11-23-51/';

%% TTLs
% get trial and odor valve transitions
temp = load(fullfile(myKsDir,'myTTLfile_1.mat'));
ch = temp.TTLs.data;
states = temp.TTLs.info.eventId;
EventTS = temp.TTLs.timestamps;

if mean(diff(EventTS))>1000
    keyboard;
    EventTS = EventTS/30000; % convert to timestamps
end
EventTS = EventTS - temp.TTLs.offset; % subtract start timestamp to align with spikes

% older recordings have channels as 0 and 1
trialCh = 1;
odorCh = 2;
if ~any(ch==2)
    trialCh = trialCh - 1;
    odorCh = odorCh - 1;
end

% Trial timestamps
TTLs.Trial(:,1) = EventTS(ch==trialCh & states,1);
while size(EventTS(ch==trialCh & ~states,1),1) > size(TTLs.Trial,1)
    TTLs.Trial(end,:) = [];
end
TTLs.Trial(:,2) = EventTS(ch==trialCh & ~states,1);
TTLs.Trial(1,:) = [];

% Odor timestamps
TTLs.Odor(:,1) = EventTS(ch==odorCh & states,1);
while size(EventTS(ch==odorCh & ~states,1),1) > size(TTLs.Odor,1)
    TTLs.Odor(end,:) = [];
end
TTLs.Odor(:,2) = EventTS(ch==odorCh & ~states,1);
if ~(TTLs.Odor(1,2) < TTLs.Trial(1,2))
    TTLs.Odor(1,:) = [];
end

%% Stimulus file
% myStimFile = fullfile(myKsDir,'250521_9_24.txt');
%
[~,ephysfile] = fileparts(myKsDir);
if nargin < 2
    % find the stim file
    serverpath = '/mnt/grid-hs/mdussauz/PhotonCerber_Stimuli_on_server';
    if ~exist(serverpath)
        keyboard;
    end

    % recording
    [~,MouseName] = fileparts(fileparts(myKsDir));
    Filename = strsplit(ephysfile,'_'); % splits date and timestamp
    FolderName = regexprep(Filename{1},'-',''); % removes dashes from date
    StimFiles = dir(fullfile(serverpath,FolderName(3:end),'*.txt'));
    if isempty(StimFiles) && strcmp(FolderName(end-1),'0')
        StimFiles = dir(fullfile(serverpath,FolderName([3:end-2 end]),'*.txt'));
    end
    if isempty(StimFiles)
        keyboard;
    else
        if size(StimFiles,1)>1
            for n = 1:size(StimFiles,1)
                Split = strsplit(StimFiles(n).name,{'_','.'});
                TimeTag(n) = str2double(Split{2})*60 + str2double(Split{3});
                numStims(n) = numel(readmatrix(fullfile(StimFiles(n).folder,StimFiles(n).name))) - 6;
            end
            TimeTagEphys = str2double(Filename{2}(1:2))*60 + str2double(Filename{2}(4:5));
            StimsEphys = size(TTLs.Trial,1);
            [deltaTime,whichFile] = min(abs(TimeTag - TimeTagEphys));
            [deltaTimeStim,whichFileStim] = min(abs(numStims - StimsEphys));
            if abs(deltaTime)>30
                keyboard;
            else
                if abs(numStims(whichFile)-StimsEphys)<=2
                    myStimFile = fullfile(StimFiles(whichFile).folder,StimFiles(whichFile).name);
                else
                    if abs(TimeTag(whichFileStim) - TimeTagEphys) < 10
                        warning(['best time matching stim file is not matching in no. of stimuli - using later file']);
                        myStimFile = fullfile(StimFiles(whichFileStim).folder,StimFiles(whichFileStim).name);
                    else
                        keyboard;
                    end
                end

            end
        else
            myStimFile = fullfile(StimFiles(1).folder,StimFiles(1).name);
        end
    end
end

% Get stimulus info from the odor machine file
StimFile = readmatrix(myStimFile);
StimSettings.timing = StimFile(1:6)'; % [whichmachine(0=16 odors) pre-stim stim post-stim iti reps]
StimFile(1:6,:) = [];
% sometimes there are trailing zeros at the end
if find(flipud(StimFile),1,'first') ~= 1
    keyboard;
    todelete = length(StimFile) - find(flipud(StimFile),1,'first') + 2;
    StimFile(todelete:end) = [];
end

% some info about number of stimuli and various reps etc
nStim = unique(StimFile);
nTypes = 2; % no. of concentrations checked
nReps = StimSettings.timing(end);
extraStim = unique(StimFile(((numel(nStim)*nReps)+1):end));
extraReps = (length(StimFile) - (numel(nStim)*nReps))/numel(extraStim);
nStim = nStim(~ismember(nStim,extraStim));

%% match the two files
if str2double(ephysfile(1:4))<2025
    % old style CID recordings my Marie on photoncerber
    switch numel(nStim)
        case 16
            SessionType = '16_Odors';
            conc_used = 10^-2;
        case 20
            SessionType = 'ConcentrationSeries';
            conc_used = [10^-4 3*10^-3 6*10^-3 10^-2];
        otherwise
            keyboard;
    end

    % make durations column
    TTLTypes = fieldnames(TTLs);
    for i = 1:numel(TTLTypes)
        TTLs.(TTLTypes{i})(:,3) = TTLs.(TTLTypes{i})(:,2) - TTLs.(TTLTypes{i})(:,1);
        TTLs.(TTLTypes{i})(TTLs.(TTLTypes{i})(:,3)<0.001,:) = [];

        splitTrials = find(abs(TTLs.(TTLTypes{i})(:,2) - circshift(TTLs.(TTLTypes{i})(:,1),-1))<0.001);
        if any(splitTrials)
            keyboard;
            disp(['merging ',num2str(numel(splitTrials)),' split ',char(TTLTypes{i}),' Trials']);
            %             Off(splitTrials,:)  = [];
            %             On(splitTrials+1,:) = [];
        end

    end

    % remove crappy trials
    trialduration = sum(StimSettings.timing(2:4)/1000);
    % first trial seems crap
    TTLs.Trial(find(TTLs.Trial(:,3)<trialduration),:) = [];
    % also next trial seems empty
    if isempty(intersect(find(TTLs.Odor(:,1)>=TTLs.Trial(1,1)),find(TTLs.Odor(:,2)<=TTLs.Trial(1,2))))
        TTLs.Odor(find(TTLs.Odor(:,1)<TTLs.Trial(2,1)),:) = [];
        TTLs.Trial(1,:) = [];
    end
    % populate TTLs.Trial with stimulus, rep, odor index
    if numel(StimFile)==(nReps*numel(nStim))
        stimlist = reshape(StimFile,[],nReps); % cols are repeats, rows are stim identities
    else
        keyboard;
    end
    count = 0;
    for rep = 1:nReps
        for stim = 1:numel(nStim)
            count = count + 1;
            if numel(conc_used)>1
                if mod(stimlist(stim,rep),5)
                    TTLs.Trial(count,4) = mod(stimlist(stim,rep),5); % which odor (divide by bank)
                else
                    TTLs.Trial(count,4) = 5;
                end
                TTLs.Trial(count,5) = conc_used(ceil(stimlist(stim,rep)/5)); % intensity
            else
                TTLs.Trial(count,4) = stimlist(stim,rep); % which odor
                TTLs.Trial(count,5) = conc_used; % intensity
            end
            TTLs.Trial(count,6) = rep; % which repeat
            % find the odor valve open and close timestamps that belong to
            % this trial
            f = intersect(find(TTLs.Odor(:,1)>=TTLs.Trial(count,1)),find(TTLs.Odor(:,2)<=TTLs.Trial(count,2)));
            if ~isempty(f)
                TTLs.Trial(count,7:9) = TTLs.Odor(f,:); % copy over the odor TTL
            else
                keyboard;
            end
        end
    end
    if ~any(TTLs.Trial(:,7)==0)
        disp('ephys trials and stimulus list matched!');
        % copy to local path
        copyfile(myStimFile,myKsDir);


        StimSettings.Odors = nStim;
        StimSettings.Reps = nReps;
        StimSettings.SessionType = SessionType;
        StimSettings.Dilutions = conc_used;
    elseif exist('todelete','var') || find(TTLs.Trial(:,7),1,'last') == length(StimFile)
        unknown = find(TTLs.Trial(:,7),1,'last') + 1;
        TTLs.Trial(unknown:end,4) = -1;

        disp('photoncerber stimulus file had missing entries!');

        % copy to local path
        copyfile(myStimFile,myKsDir);
        StimSettings.Odors = nStim;
        StimSettings.Reps = nReps;
        StimSettings.SessionType = SessionType;
        StimSettings.Dilutions = conc_used;
    end

else % new recordings that Priyanka did after adding purge etc
    if size(TTLs.Trial,1) == size(StimFile,1)
        TTLs.Trial(:,3) = TTLs.Trial(:,2) - TTLs.Trial(:,1); % Trial duration
        % add odor info to Trial TTLs
        TTLs.Trial(:,4) = StimFile; % odor identity
    else
        keyboard;
    end

    % Add odor starts and stops
    for t = 1:size(TTLs.Trial,1)
        whichrows = find(TTLs.Odor(:,1)>TTLs.Trial(t,1) & TTLs.Odor(:,2)<TTLs.Trial(t,2));
        odorTTLs = TTLs.Odor(whichrows,:)';
        TTLs.Trial(t,6+(1:numel(odorTTLs))) = odorTTLs(:);
    end
    if ~any(TTLs.Trial(:,7)==0)
        disp('ephys trials and stimulus list matched!');
        % copy to local path
        copyfile(myStimFile,myKsDir);

        StimSettings.miniOdors = nStim;
        StimSettings.miniReps = nReps;
        StimSettings.megaOdors = extraStim;
        StimSettings.megaReps = extraReps;
        StimSettings.SessionType = 'newCID';
        StimSettings.Dilutions = 0;
    end
end

if ~exist(fullfile(myKsDir,'quickprocessOdorTTLs.mat'))
    save(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings');
end

end