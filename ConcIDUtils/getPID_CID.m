function [PIDData, TTLs, StimSettings] = getPID_CID(recordNodeDir, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('whichRig','Photoncerber', @(x) ischar(x));
% params.addParameter('allUnits', 0, @(x) isnumeric(x)); % -1: non-good, 0: good (default), 1: all

% extract values from the inputParser
params.parse(varargin{:});
whichRig = 1 + strcmp(params.Results.whichRig,'Photoncerber'); % 1 = thermistor, 2 = MFS

if ~exist(fullfile(fileparts(recordNodeDir),'quickprocessPID.mat'))

    addpath(genpath('/opt/open-ephys-analysis-tools'));
    OEPSSamplingRate = 30000;

    %% from the recording directory - load the recorded PID signals
    % get aux bitVolts from the Continuous_Data.openephys
    settingsFile = fullfile(recordNodeDir, 'Continuous_Data.openephys');
    xml = readstruct(settingsFile, 'FileType', 'xml');
    whichChannel = xml.RECORDING.PROCESSOR.CHANNEL.filenameAttribute;
    auxBitVolts = xml.RECORDING.PROCESSOR.CHANNEL.bitVoltsAttribute;
    PIDFile = fullfile(recordNodeDir, whichChannel);

    % TTLs
    if ~exist(fullfile(fileparts(recordNodeDir),'myTTLfile_1.mat'))
        TTLfile = fullfile(recordNodeDir,'all_channels.events');
        [TTLs.data, TTLs.timestamps, TTLs.info] = load_open_ephys_data(TTLfile); % data has channel IDs
        NUM_HEADER_BYTES   = 1024;
        tsMap     = memmapfile(PIDFile, 'Writable', false, 'Format', 'int64', ...
            'Offset', NUM_HEADER_BYTES, 'Repeat', 1);
        startSample = tsMap.Data(1);
        startTS = double(startSample)/OEPSSamplingRate;
        TTLs.offset = startTS;
        save (fullfile(fileparts(recordNodeDir),'myTTLfile_1.mat'),'TTLs');
    end
    [TTLs,StimSettings] = makeTTLs_CID(fileparts(recordNodeDir));

    % PID - read the data 
    totalchans = 1;
    nSamples    = 1024;  % fixed to 1024 for now!
    fid{1} = fopen(PIDFile);
    % discard header information
    fseek(fid{1}, 1024, 0);

    flag = 1;
    dataOut = [];
    while 1
        samples = zeros(nSamples * 1000, totalchans, 'int16');
        for j = 1:numel(fid)
            collectSamps    = zeros(nSamples * 1000, 1, 'int16');

            rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');

            nbatches        = ceil(numel(rawData)/(nSamples+6));
            for s = 1:nbatches
                rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
            end
            samples(:,j)         = collectSamps;
        end

        if nbatches<1000
            flag = 0;
        end
        if flag==0
            samples = samples(1:s*nSamples, :);
        end

        dataOut = vertcat(dataOut,samples);


        if flag==0
            break;
        end
    end
    fclose(fid{1});

    %PID = double(dataOut)*(auxBitVolts/2) + 2.5;
    PID = double(dataOut)*(auxBitVolts);
    timestamps = (0:1:length(PID))'*(1/OEPSSamplingRate);
    timestamps(end,:) = [];

    % downsample to behavior resolution (500Hz)
    SampleRate = 500;
    PIDData(:,1) = 0:1/SampleRate:max(timestamps);
    PIDData(:,2) = interp1q(timestamps,PID,PIDData(:,1)); % thermistor

    save(fullfile(fileparts(recordNodeDir),'quickprocessPID.mat'),'TTLs','StimSettings','PIDData');
else
    load(fullfile(fileparts(recordNodeDir),'quickprocessPID.mat'),'TTLs','StimSettings','PIDData');
end

end