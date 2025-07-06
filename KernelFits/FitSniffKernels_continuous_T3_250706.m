%% script to extract the sniff aligned spiking responses
%  of a given neuron and fit the Air kernels
% settings
PSTHbinsize = 20; 
FitStyle = 'lsqcurvefit';
%KernelCondition = 3;

%% Step 1: Data Paths
myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

SavePath = '/home/priyanka/Desktop/sniffPSTHPredictions/T3_20250516';
if ~exist(SavePath,'dir')
    mkdir(SavePath);
    fileattrib(SavePath, '+w','a');
end

%% Get the sniff timestamps
%  Load Sniffs from the T3 session and make a vector with binarized sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
SniffTrace(:,1) = RespirationData(:,1); % Timestamps in 500hz base
SniffTrace(:,2) = 0;
% make a digitized trace
for i = 1:size(AllSniffs,1) % every sniff timestamp
    idx = AllSniffs(i,11:12);
    SniffTrace(idx(1):idx(2),2) = 1;
end

% figure; 
% plot(RespirationData(:,1),RespirationData(:,3));
% hold on
% plot(SniffTrace(:,1),SniffTrace(:,2));

% downsample to 20ms binsize
downsample  = PSTHbinsize/(1000/500);
TS_temp     = SniffTrace(:,1)';
TS_down     = TS_temp(1):PSTHbinsize/1000:TS_temp(end);

SniffTraceDS(:,1) = TS_down;
SniffTraceDS(:,2) = interp1(TS_temp,SniffTrace(:,2)',TS_down,'nearest');

% some extra sniff processing for later
load(fullfile(myKsDir,'SessionDetails.mat'));
startTS = 0;
for i = 1:numel(Files.Samples)
    endTS = Files.Samples(i)/30000 + startTS; % in seconds
    whichsniffs = find( (AllSniffs(:,1)>=startTS) & (AllSniffs(:,1)<endTS) );
    AllSniffs(whichsniffs,8) = i; % session phase
    startTS = endTS;
end

% group sniffs by odors
% there are no air OFF sniffs here
AllSniffs(:,4) = 1; % air is always on
[ParsedSniffs, StimulusList] = ParseSniffsByStimuli(AllSniffs, 'SortBy', 3);
Sniffs2Use{1} = ParsedSniffs{2}; % Air sniffs
Sniffs2Use{1}(:,8) = 1;

%% For fitting
InputVector = []; %zeros(6,nBins); % sniff, air, odor1, odor2, odor3, motor location
InputVector(1,:) = SniffTraceDS(1:31500,2);
InputVector(2:7,:) = 0;

%% Get the PSTHs at the same resolution
%  Load spikes
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

nBins = size(InputVector,2);
PSTHs = zeros(nUnits,nBins);
for n = 100 %:nUnits % every unit
    thisUnitSpikeTimes = SingleUnits(n).spikes; % in seconds
    thisUnitSpikeTimes = floor(thisUnitSpikeTimes*1000/PSTHbinsize) + 1; % in binIDs
    % ignore spiketimes larger than the largest sniff time
    thisUnitSpikeTimes(thisUnitSpikeTimes>nBins,:) = [];
    [C,~,ic] = unique(thisUnitSpikeTimes);
    bin_counts = accumarray(ic,1);
    if ~isempty(C)
        PSTHs(1,C) = bin_counts;
    end
end

%% set up the model
tic
telapsed = [];
for thisunit = 1 %:numel(MyUnits)
    disp(thisunit);         

    tstart = tic;

    unitid = thisunit;
    kernelLength = 700; % in ms
    StartingKernels = zeros(1,(5*(kernelLength/PSTHbinsize)) + 4); % 5 kernels + baseline + 3 locationcoeff

    Eval_max = 1e+6; Iter_max = 1e+6;
    %Fun_Tol  = 1e-8; Step_Tol = 1e-8;
    options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max); %,'TolFun',Fun_Tol,'TolX',Step_Tol);
    InputVector(8,:) = PSTHs(unitid,1:31500);
    InputVector(8,find(InputVector(7,:))) = 0; % zero out perturbation periods

    switch FitStyle
        case 'fminunc'
            %[fittedkernel{thisunit}] = fmincon(@(parameters)ssq(parameters,InputVector),StartingKernels,[],[],[],[],lb,ub,[],options);
            [fittedkernel{thisunit}] = fminunc(@(parameters)ssq(parameters,InputVector),StartingKernels,options);
        case 'fminunc_LL'
            %[fittedkernel{thisunit}] = fmincon(@(parameters)ssq(parameters,InputVector),StartingKernels,[],[],[],[],lb,ub,[],options);
            [fittedkernel{thisunit}] = fminunc(@(parameters)poisson_ll(parameters,InputVector),StartingKernels,options);            
        case 'lsqcurvefit'
            model_fit = @sniff_out; 
            lb = -100 + 0*StartingKernels; 
            ub = 100 + 0*StartingKernels;
            rectifyFR = 0; % false

            [fittedkernel{thisunit},resnorm,residual,exitflag,output] = ...
                lsqcurvefit(model_fit,StartingKernels,InputVector(1:7,:),InputVector(8,:),lb,ub,options);

    end

    telapsed(thisunit) = toc(tstart);
end

toc;

%% optimization function for lsqcurvefit
function [zdata] = sniff_out(StartingKernels,xdata)
    
    % make loal copies
    Starting_Kernels    = StartingKernels;
    x_data              = xdata;
    
    % crop out the kernels 
    %[baseline,kernels,locationcoef] = ParseSniffKernels(Starting_Kernels);
    [baseline,kernels,locationcoef] = ParseSniffKernels(Starting_Kernels,'independentcoeffs', true);
    % get psth
    [zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,x_data);
    
    %if rectifyFR
        %zdata(zdata<0) = 0;
    %end
end

%% optimization function for fminunc version 
function SSE = ssq(parameters,x)
    predictors = x(1:7,:);
    data = x(8,:);
    %datasmooth = sgolayfilt(data,1,5);
    [baseline,kernels,locationcoef] = ParseSniffKernels(parameters);
    [fitted] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,predictors);
    fitted(fitted<0) = 0;
    error1 = fitted-data;
    %error1 = fitted-datasmooth;
    SSE = sum(error1.^2);
end

%% fminunc with loglikelihood
function LLH = poisson_ll(parameters,x)
    predictors = x(1:7,:);
    data = x(8,:);
    %datasmooth = sgolayfilt(data,1,5);
    [baseline,kernels,locationcoef] = ParseSniffKernels(parameters);
    [fitted] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,predictors);
    fitted(fitted<0) = 0;
    
    LLH = -poisson_log_likelihood(data, fitted, 1);
end

function log_likelihood = poisson_log_likelihood(observed_counts, predicted_rates, bin_size)
    % Adjust the predicted rates by the bin size
    lambda = predicted_rates * bin_size * 0.99 + 0.005;
    
    % Calculate the log-likelihood for each bin
    log_likelihoods = observed_counts .* log(lambda) - lambda - gammaln(observed_counts + 1);
    
    % Sum the log-likelihoods to get the total log-likelihood
    log_likelihood = sum(log_likelihoods);
end

%%
%options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'MaxFunEvals',Eval_max,'MaxIter',Iter_max);
% 
% f = @(x)PSTHresiduals(x,InputVector,PSTHs(unitid,:));
% %KernelsOut{unitid} = fminunc(f, StartingKernels, options); 
% 
% KernelsOut{unitid} = fminsearch(f, StartingKernels); %, options); 
%toc
