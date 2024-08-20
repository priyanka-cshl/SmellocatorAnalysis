%% script to extract the sniff aligned closed-loop spiking responses
%  of a given neuron and fit the ITI, Air and Odor kernels

FitStyle = 'fminunc'; % 'fmincon' 'lsqcurvefit'
FitStyle = 'lsqcurvefit';
FitStyle = 'fminunc_LL';

%% Step 1: Get the spiking data and sniff parameters

%SessionName = 'O9_20220630_r0'; MyUnits = [];
%SessionName = 'O8_20220702_r0'; MyUnits = [];
%SessionName = 'S12_20230804_r0'; MyUnits = [];
%SessionName = 'Q4_20221109_r0'; MyUnits = [];
%SessionName = 'Q8_20221204_r0'; MyUnits = [];
%SessionName = 'Q9_20221116_r0'; MyUnits = [];
%SessionName = 'Q3_20221019_r0'; MyUnits = [];
%SessionName = 'Q4_20221109_r0'; MyUnits = [2 55 12 69 4 9 19 14 16 26 41 10]; %Q4
SessionName = 'S12_20230731_r0'; MyUnits = []; %[2 3 6 9 10 11 16 17 18 19 25 26 28 29 31 34 39 42 44 46 47 48 49 50 54 58 69 70 72 73 74 85 87 88 95 97]; % S12
%SessionName = 'Q9_20221116_r0'; MyUnits = [1 11 15 18 19 23 28 29 36 39 43 49 94]; %Q9
%SessionName = 'Q8_20221204_r0'; MyUnits = [49 54 104]; %Q8
%SessionName = 'Q3_20221019_r0'; MyUnits = [16 18]; %Q3
%SessionName = 'O3_20210927_r0'; MyUnits = [2 3 7 9 11 13 14 15 19 25 27 28 33 35 42 43 45 48 54 58 59 62 66 72] % O3

MySession = fullfile('/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/', ...
                        SessionName(1:regexp(SessionName,'_','once')-1), ...
                        [SessionName,'_processed.mat']);

SavePath = ['/home/priyanka/Desktop/sniffPSTHPredictions/', SessionName(1:regexp(SessionName,'_','once')-1)];
SavePath = fullfile(SavePath,SessionName,FitStyle);

if ~exist(SavePath,'dir')
    mkdir(SavePath);
    fileattrib(SavePath, '+w','a');
end
MatFile = fullfile(SavePath,SessionName);

[TrialAligned, TrialInfo, ...
    ReplayAligned, ReplayInfo, ...
    TuningAligned, TuningInfo, ...
    AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

if isempty(MyUnits)
    MyUnits = 1:size(AllUnits.ChannelInfo,1); % all Units 
end

%% make continous time vectors
AllSniffs = [];
for trial = 1:size(TrialAligned.Sniffs,2)
    thisTrialSniffs = TrialAligned.Sniffs{trial}(:,[7 8 2 13]); % inhalationStart inhalationEnd Odor/AirState MotorLocation
    % add the actual trial start time
    thisTrialSniffs(:,1:2) = thisTrialSniffs(:,1:2) + TrialAligned.RefTS(trial);
    AllSniffs = vertcat(AllSniffs,thisTrialSniffs);
end
% everything above is in OEPS timebase

% remove duplicates if any
if size(unique(AllSniffs(:,1:2),'rows'),1) < size(AllSniffs,1)
    keyboard;
    U = unique(AllSniffs(:,1:2),'rows');
    AllSniffs = AllSniffs(U,:);
end

% make a long vector @10ms resolution until the last sniff end + 1s
PSTHbinsize = 10; 
deltasniffs = 1; % make all sniffs the same duration

nBins = ceil(AllSniffs(end,2)*1000/PSTHbinsize);
InputVector = zeros(6,nBins); % sniff, air, odor1, odor2, odor3, motor location

for s = 1:size(AllSniffs,1) % every sniff
    ts = AllSniffs(s,1:2); % inhalation start and end in seconds 
    ts = ts*1000/PSTHbinsize; % in bins
    ts = round(ts);
    if deltasniffs
        ts(2) = ts(1);
        %ts(2) = ts(1) + (50/PSTHbinsize) - 1;
    end
    InputVector(1,ts(1):ts(2)) = 1; % inhalation period
    if AllSniffs(s,3) >= 0
        InputVector(2,ts(1):ts(2)) = 1; % Sniffs with Manifold Air ON 
        InputVector((2 + AllSniffs(s,3)),ts(1):ts(2)) = 1; % only pull up inhalations in the corresponding odor vector
    end    
    InputVector(6,ts(1):ts(2)) = AllSniffs(s,4); % keep track of odor location in each sniff
end

%% Get the PSTHs at the same resolution
PSTHs = zeros(size(AllUnits.Spikes,2),nBins);
for n = 1:size(AllUnits.Spikes,2) % every unit
    thisUnitSpikeTimes = AllUnits.Spikes{n}; % in seconds
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
    StartingKernels = zeros(1,(5*(kernelLength/PSTHbinsize)) + 2); % 5 kernels + baseline + locationcoeff

    Eval_max = 1e+6; Iter_max = 1e+6;
    %Fun_Tol  = 1e-8; Step_Tol = 1e-8;
    options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max); %,'TolFun',Fun_Tol,'TolX',Step_Tol);

    switch FitStyle
        case 'fminunc'
            InputVector(7,:) = PSTHs(unitid,:);
            %[fittedkernel{thisunit}] = fmincon(@(parameters)ssq(parameters,InputVector),StartingKernels,[],[],[],[],lb,ub,[],options);
            [fittedkernel{thisunit}] = fminunc(@(parameters)ssq(parameters,InputVector),StartingKernels,options);
        case 'fminunc_LL'
            InputVector(7,:) = PSTHs(unitid,:);
            %[fittedkernel{thisunit}] = fmincon(@(parameters)ssq(parameters,InputVector),StartingKernels,[],[],[],[],lb,ub,[],options);
            [fittedkernel{thisunit}] = fminunc(@(parameters)poisson_ll(parameters,InputVector),StartingKernels,options);
            
        case 'lsqcurvefit'
            model_fit = @sniff_out; 
            lb = -100 + 0*StartingKernels; 
            ub = 100 + 0*StartingKernels;
            %rectifyFR = 0; % false

            [fittedkernel{thisunit},resnorm,residual,exitflag,output] = ...
                lsqcurvefit(model_fit,StartingKernels,InputVector, PSTHs(unitid,:),lb,ub,options);

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
    [baseline,kernels,locationcoef] = ParseSniffKernels(Starting_Kernels);
    % get psth
    [zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,x_data);
    
%     if rectifyFR
%         zdata(zdata<0) = 0;
%     end
end

%% optimization function for fminunc version 
function SSE = ssq(parameters,x)
    predictors = x(1:6,:);
    data = x(7,:);
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
    predictors = x(1:6,:);
    data = x(7,:);
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
