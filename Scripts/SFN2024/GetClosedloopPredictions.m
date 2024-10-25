function [TS_down,PSTHs,PredictedPSTHs,WulfPSTHs] = GetClosedloopPredictions(SessionName)
% outpputs:
% TS_down: downsample timestamps for the passive period
% PSTHs: actual data (binned at the 20 ms binsize)
% PredictedPSTHs: GLM output using locally fitted kernels
% WulfPSTHs: fits for the passive part from Wulf

Paths = WhichComputer();
KernelCondition = 1; % vanilla sniff model
PSTHbinsize = 20;

%% Make the relevant input vectors to get predicted FRs 
% load the processed (wdw) session
WhereSession = fullfile(Paths.Wolf.processed,'forWDW',[SessionName,'_processed.mat']);
load(WhereSession);

%deltasniffs = 1; % make all sniffs the same duration
downsample  = PSTHbinsize/(1000/500);
TS_temp     = TracesOut.Timestamps{1};
%TS_temp     = TS_temp - TS_temp(1); % start from zero
TS_down     = TS_temp(1):PSTHbinsize/1000:TS_temp(end);

InputVector = []; %zeros(6,nBins); % sniff, air, odor1, odor2, odor3, motor location

if isfield(TracesOut,'SniffsDigitized')
    AllSniffs       = interp1(TS_temp,TracesOut.SniffsDigitized{1},TS_down,'nearest');
    OdorState       = interp1(TS_temp,TracesOut.Odor{1},TS_down,'nearest');
    ManifoldState   = interp1(TS_temp,TracesOut.Manifold{1},TS_down,'nearest');
    % parse sniffs into ITI, air, odor1, odor2, odor3
    %InputVector(1,:) = AllSniffs;

    switch KernelCondition
        case 1
            % ITI for all, Air for all odor sniffs
            InputVector(1,:) = AllSniffs; % all sniffs
            InputVector(2,:) = AllSniffs.*(OdorState>=0).*(ManifoldState==1); % air sniffs
        case 2
            % ITI for ITI only, Air for all odor sniffs
            InputVector(1,:) = AllSniffs.*(ManifoldState==0); % only ITI sniffs
            InputVector(2,:) = AllSniffs.*(OdorState>=0).*(ManifoldState==1); % air sniffs
        case 3
            % ITI for ITI only, Air for air only
            InputVector(1,:) = AllSniffs.*(ManifoldState==0); % only ITI sniffs
            InputVector(2,:) = AllSniffs.*(OdorState==0).*(ManifoldState==1); % air sniffs
    end

    InputVector(3,:) = AllSniffs.*(OdorState==1); % odor 1 sniffs
    InputVector(4,:) = AllSniffs.*(OdorState==2); % odor 2 sniffs
    InputVector(5,:) = AllSniffs.*(OdorState==3); % odor 3 sniffs
    OdorLocations = TracesOut.SniffsLocationed{1};
    OdorLocations(isnan(OdorLocations)) = 0;
    InputVector(6,:) = interp1(TS_temp,OdorLocations,TS_down); % Odor Locations

    % what to do with perturbation sniffs?
%     TrialState  = interp1(TS_temp,TracesOut.Trial{1},TS_down,'nearest');
%     InputVector(7,:) = 0; % 1 during perturbations
    InputVector(7,:) = zeros(1,size(InputVector,2)); % don't flag anything as perturbations, just predict everything

else
    disp('Curate sniffs!');
    return;    
end

%% Make PSTHs as well at the same resolution

MyUnits = 1:size(SingleUnits,2); % all Units 
nUnits = numel(MyUnits);

nBins = size(InputVector,2);
PSTHs = zeros(nUnits,nBins);
for n = 1:nUnits % every unit
    thisUnitSpikeTimes = SingleUnits(n).spikes; % in seconds
    % substract the minimum TS - time at which behavior recording start
    thisUnitSpikeTimes = thisUnitSpikeTimes - TS_down(1);
    thisUnitSpikeTimes(thisUnitSpikeTimes<0,:) = [];
    thisUnitSpikeTimes = floor(thisUnitSpikeTimes*1000/PSTHbinsize) + 1; % in binIDs
    % ignore spiketimes larger than the largest sniff time
    thisUnitSpikeTimes(thisUnitSpikeTimes>nBins,:) = [];
    [C,~,ic] = unique(thisUnitSpikeTimes);
    bin_counts = accumarray(ic,1);
    if ~isempty(C)
        PSTHs(n,C) = bin_counts;
        %PSTHs(n,C) = 1000*bin_counts/bindownto; % FR in Hz
    end
end

%% Get the predictions

for whichset = 1:6
    FitPath = fullfile(Paths.Wolf.processed,'sniffPSTHPredictions',SessionName(1:regexp(SessionName,'_','once')-1),SessionName);

    switch whichset
        case 1
            FitType = '_lsqcurvefit_1_20ms_norectify.mat';
        case 2
            FitType = '_lsqcurvefit_2_20ms_norectify.mat';
        case 3
            FitType = '_lsqcurvefit_3_20ms_norectify.mat';
        case 4
            FitType = '_lsqcurvefit_1_20ms_rectify.mat';
        case 5
            FitType = '_lsqcurvefit_2_20ms_rectify.mat';
        case 6
            FitType = '_lsqcurvefit_3_20ms_rectify.mat';
    end

    FitPath = fullfile(FitPath,[SessionName,FitType]);
    load(FitPath, 'fittedkernel','PSTHbinsize');
    %load(FitPath, 'InputVector','PSTHs','fittedkernel','PSTHbinsize');

    %% get predictions
    for i = 1:nUnits
        % get kernels
        [baseline,kernels,locationcoef] = ParseSniffKernels(fittedkernel{i},'independentcoeffs', true);
        [zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,InputVector(1:7,:));
        zdata(zdata<0) = 0;
        PredictedPSTHs{whichset}(i,:) = zdata;
    end

end

%% get fits from wulf
wulfPath = fullfile(Paths.Wolf.processed,'sniffPSTHPredictions',SessionName(1:regexp(SessionName,'_','once')-1),SessionName,[SessionName,'_wdw.mat']);
WulfPSTHs = LoadWulfGLMOutput(wulfPath,nUnits,TS_down);

end