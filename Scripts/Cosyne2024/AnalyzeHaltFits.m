%% analyzing passive replays
SessionName = 'O3_20210927_r0';

%% load the wdw processed sessions
Paths = WhichComputer();
WhereSession = fullfile(Paths.Wolf.processed,'forWDW',[SessionName,'_processed.mat']);
load(WhereSession);
binsize = mean(diff(TracesOut.Timestamps{1})); % in seconds
bufferIndices = round(0.1/binsize);
bufferSize = bufferIndices*binsize;

%% Identify the halt periods
HaltVector = TracesOut.Trial{1};
HaltVector(HaltVector>=0) = 0;
HaltVector = abs(HaltVector);
HaltStarts = find(diff(HaltVector)==1);
HaltStops = find(diff(HaltVector)==-1);
if numel(HaltStarts) == numel(HaltStops)
    HaltIdx = [HaltStarts+1 HaltStops];
else
    keyboard;
end
nHalts = size(HaltIdx,1);

%% get the predicted closed-loop PSTHs for all units
[ClosedloopFitTimestamps,ClosedloopPSTHs,PredictedClosedloopPSTHs,WulfClosedloopPSTHs] = GetClosedloopPredictions(SessionName);

%% get the predicted passive PSTHs for all units
[PassiveFitTimestamps,PassivePSTHs,PredictedPassivePSTHs,WulfPassivePSTHs] = GetPassivePredictions(SessionName);

%% Get predicted and observed FRs for every halt trial
nUnits = size(SingleUnits,2);
for whichUnit = 1:nUnits
    thisUnitSpikes = SingleUnits(whichUnit).spikes;
    
    for n = 1:nHalts % every halt trial
        ts = TracesOut.Timestamps{1}(HaltIdx(n,:)); % start and stop timestamp of halt trial
        [~,idx1] = min(abs(ClosedloopFitTimestamps-ts(1)));
        [~,idx2] = min(abs(ClosedloopFitTimestamps-ts(2)));

        myPSTH   = ClosedloopPSTHs(whichUnit,idx1:idx2); % observed
        myPSTH   = [1 numel(myPSTH) myPSTH];
        ObservedFRs{whichUnit}(n,1:numel(myPSTH)) = myPSTH;

        myPSTH   = PredictedClosedloopPSTHs{4}(whichUnit,idx1:idx2); % predicted
        myPSTH   = [1 numel(myPSTH) myPSTH];
        PredictedFRs{whichUnit}(n,1:numel(myPSTH)) = myPSTH;

        myPSTH   = WulfClosedloopPSTHs(whichUnit,idx1:idx2); % predicted
        myPSTH   = [1 numel(myPSTH) myPSTH];
        WulfFRs{whichUnit}(n,1:numel(myPSTH)) = myPSTH;
    end
end

%% some settings
dt = 0.020; % in seconds
PSTHbinsize = 20;
multfact = 1000/PSTHbinsize;
sbin = 2;

%% plot an example unit
unitIDs = cell2mat({SingleUnits.id})';
chosenUnits = 486; %[542];

for myUnit = 1:numel(chosenUnits)
    whichUnit = find(unitIDs==chosenUnits(myUnit));
    %whichUnit = find(unitIDs==chosenUnits(6));
    whichUnit = 66;
    thisUnitSpikes = SingleUnits(whichUnit).spikes;
    
    figure;

    subplot(1,3,1);
    hold on
    for n = 1:nHalts
        ts = TracesOut.Timestamps{1}(HaltIdx(n,:));
        thisHaltSpikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<=ts(2)));
        thisHaltSpikes = thisUnitSpikes(thisHaltSpikes) - ts(1);
        % plot the spikes
        PlotRaster_v2(thisHaltSpikes,n,Plot_Colors('k'),'tickwidth',1);
        % keep track of the single trial PSTHs
        %[smoothPSTH(n,:)] = MakePSTH(thisHaltSpikes',0,[0 1000*ceil(ts(2)-ts(1))],'downsample',1000/sbin,'kernelsize',100);
    end

    subplot(1,3,2);
    hold on
    for n = 1:nHalts
        myPSTH = PredictedFRs{whichUnit}(n,:);
        myPSTH = myPSTH(2+(1:myPSTH(2)));
        predictedHaltSpikes = PechePourPoisson(multfact*myPSTH,dt);
        % plot the spikes
        PlotRaster_v2(predictedHaltSpikes,n,Plot_Colors('k'),'tickwidth',1);
        % keep track of the single trial PSTHs
        %[predictedPSTH(n,:)] = MakePSTH(predictedHaltSpikes',0,[0 1000*ceil(ts(2)-ts(1))],'downsample',1000/sbin,'kernelsize',100);
    end

    subplot(1,3,3);
    hold on
    for n = 1:nHalts
        myPSTH = WulfFRs{whichUnit}(n,:);
        myPSTH = myPSTH(2+(1:myPSTH(2)));
        predictedHaltSpikes = PechePourPoisson(multfact*myPSTH,dt);
        % plot the spikes
        PlotRaster_v2(predictedHaltSpikes,n,Plot_Colors('k'),'tickwidth',1);
        % keep track of the single trial PSTHs
        %[predictedPSTH(n,:)] = MakePSTH(predictedHaltSpikes',0,[0 1000*ceil(ts(2)-ts(1))],'downsample',1000/sbin,'kernelsize',100);
    end

end

%% Quantifying fit similarity
% smoothing: sgolayfilt([],1,9);
smoothfactor = [1 9];
FitCorrelations = [];
FitResiduals    = [];
LLH             = [];
for i = 1:nUnits
    for j = 1:size(ObservedFRs{1},1)
        observed    = ObservedFRs{i}(j,:);
        observed    = observed(2+(1:observed(2)));

        predicted   = PredictedFRs{i}(j,:);
        predicted   = predicted(2+(1:predicted(2)));

        LLH(i,j) = poisson_log_likelihood(observed, predicted, 1);

        observed    = sgolayfilt(observed,smoothfactor(1),smoothfactor(2));
        fitcorr     = corrcoef(observed,predicted);
        FitCorrelations(i,j) = fitcorr(1,2);

        residual    = mean((observed - predicted).^2);
        FitResiduals(i,j) = residual;
    end

    ResidualsMean(i,1)  = mean(FitResiduals(i,12:15),'omitnan');
    ResidualsSTD(i,1)   = std(FitResiduals(i,12:15),'omitnan');
    nsamps              = 4;
    ts                  = tinv([0.025  0.975],(nsamps-1));
    ResidualsCI95(i,1)  = ts(2)*ResidualsSTD(i,1)/sqrt(nsamps);

end

%% plot an example unit
unitIDs = cell2mat({SingleUnits.id})';
chosenUnits = [921 1041 765 1218 1128 1095 979 777];
nTrials = nActiveReplays + nReplays + 1;
nTrials = 5 + 1;
FRmax = 50;
myxlims = [33 56] ; %[33 56] - 0;
figure;
nRows = numel(chosenUnits);
nCols = 4;
dt = 0.020; % in seconds
PSTHbinsize = 20;
sbin = 2;

%plotspikes = 1;

% for myUnit = 1:8
%     whichUnit = myUnit + 56;
for myUnit = 1: numel(chosenUnits)
    whichUnit = find(unitIDs==chosenUnits(myUnit));
    %whichUnit = find(unitIDs==chosenUnits(6));
    thisUnitSpikes = SingleUnits(whichUnit).spikes;
    smoothPSTH = []; fitPSTH = [];

    % plotting the template spikes
    % to get undistorted spikes
    ts = TracesOut.Timestamps{1}([TemplateIdx(1,3) TemplateIdx(end,2)]);
    % TemplateSpikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<=ts(2)));
    % TemplateSpikes = thisUnitSpikes(TemplateSpikes) - ts(1);
    % plot(TemplateSpikes,-1*ones(numel(TemplateSpikes),1),'*k');

    % to adjust spike times to account for the gaps
    TemplateSpikesAdjusted = [];
    for n = 1:size(TemplateIdx,1)
        ts = TracesOut.Timestamps{1}(TemplateIdx(n,[4 2]));
        if n == 1
            templateStart = ts(1);
        end
        thisSubtrialSpikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<=ts(2)));
        thisSubtrialSpikes = thisUnitSpikes(thisSubtrialSpikes) - templateStart;
        if n > 1
            % adjust spike times such that the extra samples are accounted for
            thisSubtrialSpikes = thisSubtrialSpikes - sum(extraindices(2:n))*binsize;
        end
        TemplateSpikesAdjusted = vertcat(TemplateSpikesAdjusted,thisSubtrialSpikes);
    end

    subplot(nRows,nCols,myUnit*nCols-3);
    % plot the trial structure
    PlotBehavior(ReplayTS,[],[],[],[],ReplayTrialVec,[],nTrials);
    set(gca,'YLim',[-0.4 nTrials+0.4],'YTick',[],'TickDir','out','XLim',[ReplayTS(1) ReplayTS(end)],'XTick', []);
    axis manual;
    hold on;
    % plot adjusted template spikes
    PlotRaster_v2(TemplateSpikesAdjusted,1,Plot_Colors('k'),'tickwidth',1);
    [smoothPSTH(1,:)] = MakePSTH(TemplateSpikesAdjusted',0,[0 1000*ceil(ts(2)-templateStart)],'downsample',1000/sbin,'kernelsize',100);

    % for the passive replays
    for n = 1:5 %nReplays % each passive replay
        ts = PassiveOut.Timestamps{1}(ReplayIdx(n,1:2));
        thisReplaySpikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<=ts(2)));
        thisReplaySpikes = thisUnitSpikes(thisReplaySpikes) - ts(1);
        % plot the spikes
        PlotRaster_v2(thisReplaySpikes,n+1,Plot_Colors('r'),'tickwidth',1);
        % keep track of the single trial PSTHs
        [smoothPSTH(n+1,:)] = MakePSTH(thisReplaySpikes',0,[0 1000*ceil(ts(2)-ts(1))],'downsample',1000/sbin,'kernelsize',100);
    end
    set(gca,'XLim',myxlims);


    % first plot the observed PSTHs
    subplot(nRows,nCols,myUnit*nCols-2);
    % plot the trial structure
    PlotBehavior(ReplayTS,[],[],[],[],ReplayTrialVec,[],FRmax);
    set(gca,'TickDir','out','XLim',[ReplayTS(1) ReplayTS(end)],'XTick', []);
    %axis manual;
    hold on;
    plot(0.001*sbin*(1:size(smoothPSTH,2)),smoothPSTH(1,:),'k','Linewidth',2);
    plot(0.001*sbin*(1:size(smoothPSTH,2)),mean(smoothPSTH(2:end-1,:),1),'r','Linewidth',2);
    set(gca,'XLim',myxlims,'YLim',[0 FRmax]);

    % predicted stuff

    if nCols>2
        % for plotting rasters
        subplot(nRows,nCols,myUnit*nCols-1);
        % plot the trial structure
        PlotBehavior(ReplayTS,[],[],[],[],ReplayTrialVec,[],nTrials);
        set(gca,'TickDir','out','XLim',[ReplayTS(1) ReplayTS(end)],'XTick', []);
        %axis manual;
        hold on;

        % for plotting PSTHs
        subplot(nRows,nCols,myUnit*nCols);
        % plot the trial structure
        PlotBehavior(ReplayTS,[],[],[],[],ReplayTrialVec,[],FRmax);
        set(gca,'TickDir','out','XLim',[ReplayTS(1) ReplayTS(end)],'XTick', []);
        %axis manual;
        hold on;

        % for the passive replays - get spike times given the predicted FR from the GLM
        % plot the closed loop
        % first predict spike times for the undistorted template stretch
        ts = TracesOut.Timestamps{1}([TemplateIdx(1,3) TemplateIdx(end,2)]); % entire set of template trials
        [~,idx1] = min(abs(ClosedloopFitTimestamps-ts(1)));
        [~,idx2] = min(abs(ClosedloopFitTimestamps-ts(2)));
        myPSTH = PredictedFRs{whichUnit}(1,:);
        myPSTH = myPSTH(2+(1:myPSTH(2)));
        predictedTemplateSpikes = PechePourPoisson(multfact*myPSTH,dt) + ts(1);

        % to adjust spike times to account for the gaps
        PredictedTemplateSpikesAdjusted = [];
        for n = 1:size(TemplateIdx,1)
            ts = TracesOut.Timestamps{1}(TemplateIdx(n,[4 2]));
            if n == 1
                templateStart = ts(1);
            end
            thisSubtrialSpikes = intersect(find(predictedTemplateSpikes>=ts(1)),find(predictedTemplateSpikes<=ts(2)));
            thisSubtrialSpikes = predictedTemplateSpikes(thisSubtrialSpikes) - templateStart;
            if n > 1
                % adjust spike times such that the extra samples are accounted for
                thisSubtrialSpikes = thisSubtrialSpikes - sum(extraindices(2:n))*binsize;
            end
            PredictedTemplateSpikesAdjusted = vertcat(PredictedTemplateSpikesAdjusted,thisSubtrialSpikes);
        end

        [fitPSTH(1,:)] = MakePSTH(PredictedTemplateSpikesAdjusted',0,[0 1000*ceil(ts(2)-templateStart)],'downsample',1000/sbin,'kernelsize',100);

        subplot(nRows,nCols,myUnit*nCols-1);
        PlotRaster_v2(PredictedTemplateSpikesAdjusted,1,Plot_Colors('k'),'tickwidth',1);
        for n = 1:5 %nReplays
            ts = PassiveOut.Timestamps{1}(ReplayIdx(n,1:2));
            % find the relevant indices in the predicts
            [~,idx1] = min(abs(PassiveFitTimestamps-ts(1)));
            [~,idx2] = min(abs(PassiveFitTimestamps-ts(2)));
            myPSTH = PredictedFRs{whichUnit}(n+nActiveReplays+1,:);
            myPSTH = myPSTH(2+(1:myPSTH(2)));
            predictedSpikes = PechePourPoisson(multfact*myPSTH,dt);
            PlotRaster_v2(predictedSpikes,n+1,Plot_Colors('r'),'tickwidth',1);
            [fitPSTH(n+1,:)] = MakePSTH(predictedSpikes',0,[0 1000*ceil(ts(2)-ts(1))],'downsample',1000/sbin,'kernelsize',100);
        end
        set(gca,'XLim',myxlims);

        subplot(nRows,nCols,myUnit*nCols);
        %area(0.001*sbin*(1:size(smoothPSTH,2)),smoothPSTH(1,:),'Facecolor',[0.5 0.5 0.5],'Linestyle','none');
        plot(0.001*sbin*(1:size(fitPSTH,2)),fitPSTH(1,:),'k','Linewidth',2);
        plot(0.001*sbin*(1:size(fitPSTH,2)),mean(fitPSTH(2:end-1,:),1),'r','Linewidth',2);
        set(gca,'XLim',myxlims,'YLim',[0 FRmax]);

    end


end
set(gcf,'Position',[1 54 1440 740]);
