function [TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves(SessionPath,MyUnits,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('binsize', 50, @(x) isnumeric(x));
params.addParameter('tuningbins', 10, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
Binsize = params.Results.binsize;
tuningbinsize = params.Results.tuningbins;
customBins = [];
switch tuningbinsize
    case 15
        TuningBinType = 'Odor';
        customBins(:,1) = (-90:15:90)' - 7.5;
        customBins(:,2) = (-90:15:90)' + 7.5;
        customBins = 125 + customBins; % make center 125
    case 10
        TuningBinType = 'Odor10';
    case 24
        TuningBinType = 'Odor24';
    otherwise
        TuningBinType = 'Odor';
end

if strcmp(computer, 'MACI64')
    datapath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior';
elseif strcmp(computer,'PCWIN64')
    datapath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

%% which Units to use
N = size(SingleUnits,2); % #units
if nargin < 2 || isempty(MyUnits)
    MyUnits = 1:N;
end

%% Get all spikes, all units aligned to trials - both for closed loop and replays
% hack: add a fake field TargetEntry to TrialInfo so that
% TrialAlignedSpikes doesn't crash
TrialInfo.TargetEntry = NaN*TrialInfo.Odor;

[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
                            size(TrialInfo.TrialID,2),TrialInfo,MySession);

% Events: OdorStart, Reward, TrialOFF, First TZ entry
BehaviorBinsize = Binsize/(1000/SampleRate);

%% Assembling tuning curves for non-perturbed, closed loop trials
for whichodor = 1:3
    % 1. Pick the right trials - no perturbation, and a given odor
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
        find(TrialInfo.Odor==whichodor));
    XVar = []; YVar = []; TossVar = [];
    for whichUnit = 1:numel(MyUnits) % for every Unit
        counts = [1 0];
        thisUnitSpikes = AlignedSpikes(:,MyUnits(whichUnit));
        for x = 1:size(whichTrials,1) % every trial
            % only keep the PSTH and behavioral variables from odor start until Trial OFF
            t1 = round((startoffset + Events(whichTrials(x,1),1))*SampleRate); % odor start
            t1 = max(1,t1);
            t2 = round(Events(whichTrials(x,1),3)*SampleRate); % trial OFF
            if mod(numel(t1:t2),BehaviorBinsize)
                t2 = t1 + BehaviorBinsize*floor(numel(t1:t2)/BehaviorBinsize) - 1;
            end
            newSamples = numel(t1:t2)/BehaviorBinsize;
            counts(2) = counts(2) + newSamples; % one extra bin for pre-odor start spikes
            if whichUnit == 1
                % get behavioral variable
                myMotor = mean(reshape(Traces.Motor{whichTrials(x,1)}(t1:t2),BehaviorBinsize,[]));
                XVar(counts(1):counts(2),1) = myMotor';
                TossVar(counts(1):counts(2),1) = round(rand(1));
            end
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1}; % aligned to trial ON
            thisTrialSpikes = thisTrialSpikes - Events(whichTrials(x,1)); % subtract odor start - odor start becomes zero
            % convert to ms
            thisTrialSpikes = ceil(thisTrialSpikes*1000/Binsize); % at binned resolution
            thisTrialSpikes(thisTrialSpikes<1) = [];
            [C,~,ic] = unique(thisTrialSpikes);
            bin_counts = accumarray(ic,1);
            myRaster = [];
            if ~isempty(C)
                myRaster(C) = bin_counts;
            else
                myRaster = zeros(1,newSamples);
            end
            if numel(myRaster)<diff(counts)+1
                myRaster = [myRaster, zeros(1,(diff(counts) + 1 - numel(myRaster)))];
            end
            YVar(counts(1):counts(2),whichUnit) = myRaster(1:(1+diff(counts)));
            
            counts(1) = counts(2)+1;
        end
    end
    [myCurve, XBins] = SmellocatorTuning(TuningBinType,125-XVar,1000*(YVar/Binsize),customBins);
    TuningCurve.ClosedLoopFull(:,:,:,whichodor) = squeeze(myCurve(:,[1 3],:,:));
    [myCurve, ~] = SmellocatorTuning(TuningBinType,125-XVar(find(TossVar),:),1000*(YVar(find(TossVar),:)/Binsize),customBins);
    TuningCurve.ClosedLoopHalf1(:,:,:,whichodor) = squeeze(myCurve(:,[1 3],:,:));
    [myCurve, ~] = SmellocatorTuning(TuningBinType,125-XVar(find(~TossVar),:),1000*(YVar(find(~TossVar),:)/Binsize),customBins);
    TuningCurve.ClosedLoopHalf2(:,:,:,whichodor) = squeeze(myCurve(:,[1 3],:,:));
end


PairedCorrs = [];
PairedResiduals = [];
ControlCorrs = [];
ControlResiduals = [];
%% tuning curves for replay trials                                          
if any(strcmp(TrialInfo.Perturbation,'OL-Replay'))
    [ReplayAlignedSpikes, ReplayEvents, ReplayInfo] = ...
        ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,TrialInfo,Events);
    
TemplateTrials = find(strcmp(TrialInfo.Perturbation,'OL-Template'));

    for whichodor = 1:3
        % 1. Pick the right trials - no perturbation, and a given odor
        whichTrials = find(ReplayInfo.Odor==whichodor);
        XVar = []; YVar = [];
        countsactive = 0;
        for whichUnit = 1:numel(MyUnits) % for every Unit
            counts = [1 0];
            thisUnitSpikes = ReplayAlignedSpikes(:,MyUnits(whichUnit));
            for x = 1:size(whichTrials,1) % every trial
                % only keep the PSTH and behavioral variables from odor start until Trial OFF
                t1 = round((startoffset + ReplayEvents(whichTrials(x,1),1))*SampleRate);
                t2 = round(ReplayEvents(whichTrials(x,1),3)*SampleRate);
                if mod(numel(t1:t2),Binsize/2)
                    t2 = t1 + (Binsize/2)*floor(numel(t1:t2)/(Binsize/2)) - 1;
                end
                newSamples = numel(t1:t2)/(Binsize/2);
                counts(2) = counts(2) + newSamples;
                if whichUnit == 1
                    whichTemplateTrial = round(100*rem(ReplayInfo.TrialID(whichTrials(x,1)),1));
                    if whichTemplateTrial<0 && countsactive == 0
                        countsactive = counts(1) - 1;
                    end
                    % get behavioral variable
                    myMotor = mean(reshape(Traces.Motor{TemplateTrials(abs(whichTemplateTrial))}(t1:t2),Binsize/2,[]));
                    XVar(counts(1):counts(2),1) = myMotor';                    
                end
                thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1}; % aligned to trial ON
                thisTrialSpikes = thisTrialSpikes - ReplayEvents(whichTrials(x,1)); % subtract odor start
                % convert to ms
                thisTrialSpikes = ceil(thisTrialSpikes*1000/Binsize); % at binned resolution
                thisTrialSpikes(thisTrialSpikes<1) = [];
                [C,~,ic] = unique(thisTrialSpikes);
                bin_counts = accumarray(ic,1);
                myRaster = [];
                if ~isempty(C)
                    myRaster(C) = bin_counts;
                else
                    myRaster = zeros(1,newSamples);
                end
                if numel(myRaster)<diff(counts)+1
                    myRaster = [myRaster, zeros(1,(diff(counts) + 1 - numel(myRaster)))];
                end
                YVar(counts(1):counts(2),whichUnit) = myRaster(1:(1+diff(counts)));
            
                counts(1) = counts(2)+1;
            end
        end

        [myCurve, ~] = SmellocatorTuning(TuningBinType,125-XVar(1:countsactive,:),1000*(YVar(1:countsactive,:)/Binsize),customBins);
        TuningCurve.OpenLoop(:,:,:,whichodor) = squeeze(myCurve(:,[1 3],:,:));
        [myCurve, ~] = SmellocatorTuning(TuningBinType,125-XVar((countsactive+1):end,:),1000*(YVar((countsactive+1):end,:)/Binsize),customBins);
        TuningCurve.Passive(:,:,:,whichodor) = squeeze(myCurve(:,[1 3],:,:));
    end
    
    %% if there are replay trials - compare tuning curves across these conditions
    [PairedCorrs,PairedResiduals,ControlCorrs,ControlResiduals] = CompareTuningCurves(TuningCurve);
end

if ~isempty(customBins)
    XBins = XBins - 125;
end

end
