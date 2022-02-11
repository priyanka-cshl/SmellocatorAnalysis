MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat';
MySession = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/O3/O3_20211005_r0_processed.mat';

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

% hack: add a fake field TargetEntry to TrialInfo so that
% TrialAlignedSpikes doesn't crash
TrialInfo.TargetEntry = NaN*TrialInfo.Odor;

%% Get all spikes, all units aligned to trials - both for closed loop and replays
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
                            size(TrialInfo.TrialID,2),TrialInfo);

if any(strcmp(TrialInfo.Perturbation,'OL-Replay'))
    [ReplayAlignedSpikes, ReplayEvents, ReplayInfo] = ...
                    ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
                        ReplayTTLs,TrialInfo,Events);
end

N = size(SingleUnits,2); % #units
MyUnits = (1:N);
MyUnits = [58 35 34 55 21];

%% Assembling tuning curves for closed loop trials
for whichodor = 1:3
    % 1. Pick the right trials - no perturbation, and a given odor
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
        find(TrialInfo.Odor==whichodor));
    XVar = []; YVar = [];
    for whichUnit = 1:numel(MyUnits) % for every Unit
        counts = [1 0];
        thisUnitSpikes = AlignedSpikes(:,MyUnits(whichUnit));
        for x = 1:size(whichTrials,1) % every trial
            % only keep the PSTH and behavioral variables from odor start until Trial OFF
            t1 = round((startoffset + Events(whichTrials(x,1),1))*SampleRate);
            t2 = round(Events(whichTrials(x,1),3)*SampleRate);
            counts(2) = counts(2) + numel(t1:t2);
            if whichUnit == 1
                % get behavioral variable
                XVar(counts(1):counts(2),1) = Traces.Motor{whichTrials(x,1)}(t1:t2);
            end
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1}; % aligned to trial ON
            thisTrialSpikes(thisTrialSpikes<-startoffset) = []; % ignore any spikes 1 sec preceding Trial ON
            % add startoffset to make all spike times positive
            thisTrialSpikes = thisTrialSpikes + startoffset;
            % make a PSTH
            thisTrialPSTH = MakePSTH(thisTrialSpikes, 0, [], 'kernelsize',100,'downsample', SampleRate);
            % PSTH can be empty or smaller than the trialwindow of interest
            if isempty(thisTrialPSTH)
                YVar(counts(1):counts(2),whichUnit) = 0;
            elseif numel(thisTrialPSTH)<t2
                thisTrialPSTH(end+1:t2,1) = 0;
                YVar(counts(1):counts(2),whichUnit) = thisTrialPSTH(t1:t2);
            else
                YVar(counts(1):counts(2),whichUnit) = thisTrialPSTH(t1:t2);
            end
            counts(1) = counts(2)+1;
        end
    end
%     % Hack to try many different delays
%     for delays = 10:10:100
%         XVar(:,end+1) = circshift(XVar(:,1),-delays);
%     end
    % add one column for randomized locations
    XVar(:,2) = XVar(randperm(size(XVar,1)),1);
    [Curve_CL{whichodor}, XBins] = SmellocatorTuning('Odor',125-XVar,YVar);
end

%% OpenLoop and Passive Replays
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
            counts(2) = counts(2) + numel(t1:t2);
            if whichUnit == 1
                whichTemplateTrial = round(100*rem(ReplayInfo.TrialID(whichTrials(x,1)),1));
                if whichTemplateTrial<0 & countsactive == 0
                    countsactive = counts(1) - 1;
                end
                % get behavioral variable
                XVar(counts(1):counts(2),1) = Traces.Motor{abs(whichTemplateTrial)}(t1:t2);
                
            end
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1}; % aligned to trial ON
            thisTrialSpikes(thisTrialSpikes<-startoffset) = []; % ignore any spikes 1 sec preceding Trial ON
            % add startoffset to make all spike times positive
            thisTrialSpikes = thisTrialSpikes + startoffset;
            % make a PSTH
            thisTrialPSTH = MakePSTH(thisTrialSpikes', 0, [], 'kernelsize',100,'downsample', SampleRate);
            % PSTH can be empty or smaller than the trialwindow of interest
            if isempty(thisTrialPSTH)
                YVar(counts(1):counts(2),whichUnit) = 0;
            elseif numel(thisTrialPSTH)<t2
                thisTrialPSTH(end+1:t2,1) = 0;
                YVar(counts(1):counts(2),whichUnit) = thisTrialPSTH(t1:t2);
            else
                YVar(counts(1):counts(2),whichUnit) = thisTrialPSTH(t1:t2);
            end
            counts(1) = counts(2)+1;
        end
    end
%     % Hack to try many different delays
%     for delays = 10:10:100
%         XVar(:,end+1) = circshift(XVar(:,1),-delays);
%     end
    foo = XVar(1:countsactive,1);
    XVar(1:countsactive,2) = foo(randperm(numel(foo)));
    foo = XVar((countsactive+1):end,1);
    XVar((countsactive+1):end,2) = foo(randperm(numel(foo)));
    [Curve_OL{whichodor}, XBins] = SmellocatorTuning('Odor',125-XVar(1:countsactive,:),YVar(1:countsactive,:));
    [Curve_PR{whichodor}, XBins] = SmellocatorTuning('Odor',125-XVar((countsactive+1):end,:),YVar((countsactive+1):end,:));
end

%%
figure;
map = brewermap(15,'*GnBu');
for i = 1:numel(MyUnits)
    for odor = 1:3
        subplot(3,5,i+5*(odor-1)); hold on
        line([0 0],[0 20],'LineStyle',':','Color','k');
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,1,i)',Curve_CL{odor}(:,3,1,i)',Plot_Colors('k'),{},0.5);
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_CL{odor}(:,1,2,i)',Curve_CL{odor}(:,3,2,i)',Plot_Colors('k'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,1,i)',Curve_OL{odor}(:,3,1,i)',Plot_Colors('r'),{},0.5);
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_OL{odor}(:,1,2,i)',Curve_OL{odor}(:,3,2,i)',Plot_Colors('r'));
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,1,i)',Curve_PR{odor}(:,3,1,i)',Plot_Colors('t'),{},0.5);
        MyShadedErrorBar(mean(XBins,2)'-125,Curve_PR{odor}(:,1,2,i)',Curve_PR{odor}(:,3,2,i)',Plot_Colors('t'));
%         for delays = 2:size(XVar,2)
%           MyShadedErrorBar(mean(XBins,2)',Curve_CL{odor}(:,1,delays,i)',Curve_CL{odor}(:,3,delays,i)',map(delays,:));
% %         MyShadedErrorBar(mean(XBins,2)',Curve_OL(:,1,odor,i)',Curve_OL(:,3,odor,i)',Plot_Colors('r'));
% %         if odor < 4
% %             MyShadedErrorBar(mean(XBins,2)',Curve_PR(:,1,odor,i)',Curve_PR(:,3,odor,i)',Plot_Colors('t'));
% %         end
%         end
    end
end