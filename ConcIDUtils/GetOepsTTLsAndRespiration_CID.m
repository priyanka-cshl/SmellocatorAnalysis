function [TTLs,ReplayTTLs,EphysTuningTrials,SniffData] = GetOepsTTLsAndRespiration_CID(myOEPSDir, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('ADC', true, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
GetAux = params.Results.ADC;

%% defaults
Paths = WhichComputer();
addpath(genpath(fullfile(Paths.Code,'open-ephys-analysis-tools')));
addpath(genpath(fullfile(Paths.Code,'afterphy')));

%% Get Trial Timestamps from the OpenEphys Events file
filename = fullfile(myOEPSDir,'Record*','all_channels.events');

% hack to avoid going through empty files
temp = dir(filename);
if ~isempty(temp) % old GUI
    if ~temp.bytes
        TTLs = [];
        EphysTuningTrials = [];
        SniffData = [];
        ReplayTTLs = [];
        disp('empty events file');
        return
    else
        filename = fullfile(temp(1).folder,temp(1).name);
    end

    [data, timestamps, info] = load_open_ephys_data(filename); % data has channel IDs

    % adjust for clock offset between open ephys and kilosort
    [offset] = AdjustClockOffset(fileparts(filename));
else % new GUI
%     session = Session(fileparts(myOEPSDir));
%     Events = session.recordNodes{1}.recordings{1}.ttlEvents('Acquisition_Board-100.Rhythm Data');
%     data            = Events.channel;
%     timestamps      = Events.timestamp;
%     info.eventId    = Events.state;
% 
%     % adjust for clock offset between open ephys and kilosort
%     offset = session.recordNodes{1}.recordings{1}.continuous('Acquisition_Board-100.Rhythm Data').timestamps(1);
end


timestamps = timestamps - offset;

%% Get various events
TTLTypes = unique(data);
Tags = {'Trial', 'Odor'}; 
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
    SniffData = [];
    disp('no valid oeps file found');
    return
end

%% get the stimulus list (either available locally, or from the server)
[~,stimfile_token] = fileparts(myOEPSDir);
stimfile_token = [regexprep(stimfile_token(3:13),'-',''),'*.txt'];

if ~isempty(dir(fullfile(myOEPSDir,stimfile_token)))
    stimfile = dir(fullfile(myOEPSDir,stimfile_token));
elseif ~isempty(dir(fullfile(Paths.CID.PhotonCerberStimFiles,stimfile_token(1:6),stimfile_token)))
    stimfile = dir(fullfile(Paths.CID.PhotonCerberStimFiles,stimfile_token(1:6),stimfile_token));
    % copy to local path
    copyfile(fullfile(stimfile.folder,stimfile.name),myOEPSDir);
    stimfile = dir(fullfile(myOEPSDir,stimfile_token));
else
    stimfile = [];
    disp('no stim file located');
end

if ~isempty(stimfile)
    % read the sequence of stimuli and settings - add the relevant info to
    % the TTLs struct directly
    fid = fopen(fullfile(stimfile.folder,stimfile.name)); 
    stimlist = fscanf(fid,'%d'); 
    fclose(fid);
    stimsettings = stimlist(1:6);
    stimlist(1:6,:) = []; % ? pre stim post iti #reps
    trialduration = sum(stimsettings(2:4)/1000);
    stimlist = reshape(stimlist,[],stimsettings(6)); % cols are repeats, rows are stim identities
    conc_used = [10^-4 3*10^-3 6*10^-3 10^-2];

    % first trial seems crap
    TTLs.Trial(find(TTLs.Trial(:,3)<trialduration),:) = []; 
    % also next trial seems empty
    if isempty(intersect(find(TTLs.Odor(:,1)>=TTLs.Trial(1,1)),find(TTLs.Odor(:,2)<=TTLs.Trial(1,2))))
        TTLs.Odor(find(TTLs.Odor(:,1)<TTLs.Trial(2,1)),:) = [];
        TTLs.Trial(1,:) = [];
    end
    count = 0;
    for rep = 1:size(stimlist,2)
        for stim = 1:size(stimlist,1)
            count = count + 1;
            TTLs.Trial(count,4) = mod(stimlist(stim,rep),5);
            TTLs.Trial(count,5) = conc_used(ceil(stimlist(stim,rep)/5));
            TTLs.Trial(count,6) = rep;
            % find the odor valve open and close timestamps that belong to
            % this trial
            f = intersect(find(TTLs.Odor(:,1)>=TTLs.Trial(count,1)),find(TTLs.Odor(:,2)<=TTLs.Trial(count,2)));
            if ~isempty(f)
                TTLs.Trial(count,7:9) = TTLs.Odor(f,:);
            else
                keyboard;
            end
        end
    end
end

if ~any(TTLs.Trial(:,7)==0)
    disp('ephys trials and stimulus list matched!');
end

SniffData = [];
if GetAux
    %% Get analog/digital AuxData from Oeps files - for extracting sniffing
    % for reading the aux channels
    % check if there are ADC files
    myFiles = dir(fullfile(fileparts(filename),'*_ADC*.continuous'));
    if isempty(myFiles)
        myFiles = dir(fullfile(fileparts(filename),'*.continuous'));
        filetags = cellfun(@(x) str2double(x(regexpi(x,'_')+1:regexpi(x,'.continuous')-1)), {myFiles.name})';
        filetags(:,2) = 1:size(filetags,1);
        filetags = sortrows(filetags,1);
        
        filefound = 0; fromlast = 7;
        while ~filefound

            % lets assume last 8 channels are ADCs
            ADCfilename = fullfile(myFiles(1).folder,myFiles(filetags(end-fromlast,2),:).name);
            [AuxVals, Auxtimestamps, ~] = load_open_ephys_data(ADCfilename); % data has channel IDs
            Auxtimestamps = Auxtimestamps - offset;
            % downsample to behavior resolution (500 hz)
            SniffData(:,1) = 0:1/500:max(Auxtimestamps);
            SniffData(:,2) = interp1q(Auxtimestamps,AuxVals,SniffData(:,1));
            figure;
            plot(SniffData(:,1),SniffData(:,2));
            set(gca,'XLim',[5 15]);

            answer = questdlg('Does the respiration data look ok?', ...
            	'Respiration check', ...
            	'Yes','No','Yes');
            % Handle response
            switch answer
                case 'Yes'
                    filefound = 1;
                    hold on
                case 'No'
                    close(gcf);
                    if fromlast
                        fromlast = fromlast - 1;
                    else
                        disp('no respiration file found');
                        keyboard;
                    end
            end
        end
    end

    % if respiration file was found - filter and detect peaks and valleys?
    [SniffsFiltered] = FilterThermistor(SniffData(:,2));
    plot(SniffData(:,1),SniffsFiltered,'r');

end

end
