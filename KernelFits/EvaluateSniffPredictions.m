%function [] = EvaluateSniffPredictions()

% Load the predicted kernels 
FitPath = '/home/priyanka/Desktop/sniffPSTHPredictions/S12/S12_20230731_r0/S12_20230731_r0_lsqcurvefit_20ms.mat';
load(FitPath);
% loads ('SessionName','MyUnits','InputVector','PSTHs','fittedkernel','PSTHbinsize')

%%
StimVector = InputVector(1:7,:);
% include predictions in the perturbation period as well
PerturbationBins = find(StimVector(:,7));
StimVector(PerturbationBins,7) = 0;

nUnits = size(MyUnits,2);
for i = 1:nUnits
    whichunit = MyUnits(i);
    % get kernels etc
    [baseline,kernels,locationcoef] = ParseSniffKernels(fittedkernel{whichunit});

    [PredictedPSTH] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,StimVector);
    % this prediction has -ve rates as well

    Fitted.Baseline(i) = baseline;
    Fitted.Kernels(i) = {kernels};
    Fitted.LocationScaling = {locationcoef};
    Fitted.PSTH(i,:) = PredictedPSTH;
end

%%
% Get sniff starts and club by stimulus type
SniffVector = InputVector(1:6,:);
SniffVector(:,PerturbationBins) = 0; % ignore any perturbation sniffs
AllSniffs = find(diff(SniffVector(1,:)==1));
%Sniffs{1} = AllSniffs(find(~ismember(AllSniffs,find(diff(SniffVector(2,:)==1)))));
Sniffs{1} = find(diff(SniffVector(1,:)==1));
Sniffs{2} = find(diff(SniffVector(3,:)==1));
Sniffs{3} = find(diff(SniffVector(4,:)==1));
Sniffs{4} = find(diff(SniffVector(5,:)==1));

%%
PSTHWindow = 0:1:(700/PSTHbinsize);
for i = 1:nUnits
    [AvgFitPSTH] = MakeSniffAlignedPSTH(Sniffs,Fitted.PSTH(i,:),PSTHWindow);
    
end

function [AvgPSTH,MyPSTH] = MakeSniffAlignedPSTH(Sniffs,FitPSTHs,PSTHWindow)
    if nargin < 3
        PSTHWindow = [0:1:70];
    end

    % sniff aligned psths
    for stimtype = 1:4
        PSTHtemp = [];
        AllSniffs = Sniffs{stimtype}(1:end);
        for n = 1:numel(AllSniffs)
            if (AllSniffs(n) + PSTHWindow(end)) <= numel(FitPSTHs)
                PSTHtemp(n,:) = FitPSTHs(1,AllSniffs(n)+PSTHWindow);
            end
        end
        AvgPSTH{stimtype} = mean(PSTHtemp,1);
        MyPSTH{stimtype} = PSTHtemp;
    end
end
%end