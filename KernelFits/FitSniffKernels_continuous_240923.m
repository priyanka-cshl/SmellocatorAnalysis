%% script to extract the sniff aligned closed-loop spiking responses
%  of a given neuron and fit the ITI, Air and Odor kernels

% FitStyle = 'fminunc'; % 'fmincon' 'lsqcurvefit'
% FitStyle = 'fminunc_LL';
FitStyle = 'lsqcurvefit';
KernelCondition = 3;

%% Step 1: Data Paths
SessionName = 'S12_20230731_r0'; 
MyUnits = []; % just a hack for processing all units

% MySession = fullfile('/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/', ...
%                         SessionName(1:regexp(SessionName,'_','once')-1), ...
%                         [SessionName,'_processed.mat']);

MySession = fullfile('/home/priyanka/Desktop/forWDW',[SessionName,'_processed.mat']);
SavePath = ['/home/priyanka/Desktop/sniffPSTHPredictions/', SessionName(1:regexp(SessionName,'_','once')-1)];
SavePath = fullfile(SavePath,SessionName);
if ~exist(SavePath,'dir')
    mkdir(SavePath);
    fileattrib(SavePath, '+w','a');
end
MatFile = fullfile(SavePath,[SessionName,'_',FitStyle,'_',num2str(KernelCondition),'_20ms_norectify.mat']);

%% Loading the actual data
% use the same rules as for preprocessing for wolf
%[TracesOut, SingleUnits] = LoadLeverTaskSession_wdw(MySession,0);
load(MySession);

% settings
PSTHbinsize = 20; 
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
    TrialState  = interp1(TS_temp,TracesOut.Trial{1},TS_down,'nearest');
    InputVector(7,:) = (TrialState<0); % 1 during perturbations

else
    disp('Curate sniffs!');
    return;
end

%% Get the PSTHs at the same resolution

if isempty(MyUnits)
    MyUnits = 1:size(SingleUnits,2); % all Units 
end
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

%% set up the model
tic
telapsed = [];
for thisunit = 1:numel(MyUnits)
    disp(thisunit);         

    tstart = tic;

    unitid = MyUnits(thisunit);
    kernelLength = 700; % in ms
    %StartingKernels = zeros(1,(5*(kernelLength/PSTHbinsize)) + 2); % 5 kernels + baseline + locationcoeff
    StartingKernels = zeros(1,(5*(kernelLength/PSTHbinsize)) + 4); % 5 kernels + baseline + 3 locationcoeff

    Eval_max = 1e+6; Iter_max = 1e+6;
    %Fun_Tol  = 1e-8; Step_Tol = 1e-8;
    options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max); %,'TolFun',Fun_Tol,'TolX',Step_Tol);
    InputVector(8,:) = PSTHs(unitid,:);
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

%%
save(MatFile,'SessionName','MyUnits','InputVector','PSTHs','fittedkernel','PSTHbinsize');

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
        zdata(zdata<0) = 0;
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
