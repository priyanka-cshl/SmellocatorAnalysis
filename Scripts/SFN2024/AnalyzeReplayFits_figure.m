%% analyzing passive replays
SessionName = 'O3_20211005_r0';

%% load the wdw processed sessions
Paths = WhichComputer();
WhereSession = fullfile(Paths.Wolf.processed,'forWDW',[SessionName,'_processed.mat']);
load(WhereSession);
binsize = mean(diff(TracesOut.Timestamps{1})); % in seconds
bufferIndices = round(0.1/binsize);
bufferSize = bufferIndices*binsize;

%% identify the passive replay stretches
TrialVector = PassiveOut.Replay{1};
TrialVector(TrialVector>-10) = 0;
replayedTrials = unique(TrialVector(find(TrialVector)));
nTemplates = numel(find(diff(replayedTrials)~=1)) + 1;

if nTemplates == 1
    subtrials = sort(abs(replayedTrials));
    % find the stretches of replaym
    ReplayIdx(:,1) = find(diff(PassiveOut.Replay{1})==-subtrials(1)) + 1;
    ReplayIdx(:,2) = find(diff(PassiveOut.Replay{1})==subtrials(end));

    %     figure;
    %     hold on;
    %     for nReplays = 1:size(ReplayIdx,1)
    % %         % plotting Motor
    % %         plot(PassiveOut.Motor{1}(ReplayIdx(nReplays,1):ReplayIdx(nReplays,2)));
    %
    %         % plotting sniffs
    %         subplot(size(ReplayIdx,1),1,nReplays);
    %         plot(PassiveOut.SniffsFiltered{1}(ReplayIdx(nReplays,1):ReplayIdx(nReplays,2)));
    %     end

    % for every replay: find the subtrial start and stop indices
    nReplays = size(ReplayIdx,1);
    for whichReplay = 1:nReplays
        myIndices = ReplayIdx(whichReplay,1):ReplayIdx(whichReplay,2);
        thisReplaystretch = PassiveOut.Replay{1}(myIndices);
        subtrialIdx = [];
        subtrialIdx(:,1) = vertcat(0, find(diff(thisReplaystretch)<0)) +  1;
        subtrialIdx(:,2) = vertcat(find(diff(thisReplaystretch)>0), numel(thisReplaystretch));
        if whichReplay == 1
            ReplayOdorSeq   = PassiveOut.Odor{1}(subtrialIdx(:,1)+ReplayIdx(whichReplay,1));
            ReplayTS        = PassiveOut.Timestamps{1}((myIndices(1)-bufferIndices):(myIndices(end)+bufferIndices));
            ReplayTS        = ReplayTS - PassiveOut.Timestamps{1}(myIndices(1));
            ReplayTrialVec  = PassiveOut.Odor{1}((myIndices(1)-bufferIndices):(myIndices(end)+bufferIndices));
        end
        ReplayIdx(whichReplay,2+(1:numel(subtrialIdx))) = reshape(subtrialIdx',numel(subtrialIdx),[])';
    end

end

SubtrialDurations = mode(ReplayIdx(:,4:2:end) - ReplayIdx(:,3:2:end))'; % can be used for template checking
SubTrialGaps = -mode([zeros(nReplays,1) ReplayIdx(:,4:2:end-2)] - ReplayIdx(:,4:2:end))'; % this will be used to figure out how many samples were dropped at each trial end from the CL template

%% Compare to template and find the relevant matching indices?
TrialIdxCL(:,1) = find(diff([0; abs(TracesOut.Trial{1})])>0); % trial start (not odor start)
TrialIdxCL(:,2) = find(diff([abs(TracesOut.Trial{1}); 0])<0); % trial end

TemplateIdx = TrialIdxCL(subtrials,:); % trial start and end of the trials used as template
TemplateOdorSeq = TracesOut.Odor{1}(TemplateIdx(:,1)+1); % odor sequence
if ~isequal(TemplateOdorSeq,ReplayOdorSeq) % does it match with the replays
    disp('replay mismatch');
    keyboard;
end

% get the odorstart time instead of trial start time to compare with the replay subtrials
for n = 1:size(TemplateIdx,1)
    TemplateIdx(n,3) = find(TracesOut.Odor{1}(1:TemplateIdx(n,2))~=TemplateOdorSeq(n),1,'last') + 1;
end

% get trial-to-trial duration to compute extra template samples
TemplateTrialGaps = diff([TemplateIdx(1,1); TemplateIdx(:,2)]);
extraindices = TemplateTrialGaps - SubTrialGaps;
% which indices post every trial off did the replay subtrial actually begin
TemplateIdx(2:end,4) = TemplateIdx(1:end-1,2) + extraindices(2:end) + 1;
TemplateIdx(1,4) = TemplateIdx(1,1) + extraindices(1) + 1;

% concatenate all valid indices
TemplateIndices = [];
for n = 1:size(TemplateIdx,1)
    TemplateIndices = horzcat(TemplateIndices, TemplateIdx(n,4):TemplateIdx(n,2));
end

% check the match (use duration)
if any(abs(ReplayIdx(:,2)-ReplayIdx(:,1)-numel(TemplateIndices))>2)
    disp('replay mismatch');
    keyboard;
end

total_indices = min([numel(TemplateIndices); ReplayIdx(:,2)-ReplayIdx(:,1)]);

%% Identify the active replay stretches
TVec = TracesOut.Trial{1};
TVec(TVec>0) = 0;
TVec = abs(TVec);
TrialIdxOL(:,1) = find(diff([0; TVec])>0); % trial start (not odor start)
TrialIdxOL(:,2) = find(diff([TVec; 0])<0); % trial end

nActiveReplays = size(TrialIdxOL,1);
%figure; hold on
% for m = 1:nActiveReplays
%     % adjust start and ends to begin with the odor start and odor end
%     plot(TracesOut.Motor{1}(TrialIdxOL(m,1):TrialIdxOL(m,2)));
% end

%% get the predicted closed-loop PSTHs for all units
[ClosedloopFitTimestamps,ClosedloopPSTHs,PredictedClosedloopPSTHs,WulfClosedloopPSTHs] = GetClosedloopPredictions(SessionName);

%% get the predicted passive PSTHs for all units
[PassiveFitTimestamps,PassivePSTHs,PredictedPassivePSTHs,WulfPassivePSTHs] = GetPassivePredictions(SessionName);

%% Get predicted and observed FRs for every replay and the template
nUnits = size(SingleUnits,2);
for whichUnit = 1:nUnits
    thisUnitSpikes = SingleUnits(whichUnit).spikes;

    % for the template
    ts = TracesOut.Timestamps{1}([TemplateIdx(1,3) TemplateIdx(end,2)]); % start and stop timestamp of the entire set of template trials
    [~,idx1] = min(abs(ClosedloopFitTimestamps-ts(1)));
    [~,idx2] = min(abs(ClosedloopFitTimestamps-ts(2)));

    myPSTH   = ClosedloopPSTHs(whichUnit,idx1:idx2); % observed
    myPSTH   = [1 numel(myPSTH) myPSTH];
    ObservedFRs{whichUnit}(1,1:numel(myPSTH)) = myPSTH;

    myPSTH   = PredictedClosedloopPSTHs{4}(whichUnit,idx1:idx2); % predicted
    myPSTH   = [1 numel(myPSTH) myPSTH];
    PredictedFRs{whichUnit}(1,1:numel(myPSTH)) = myPSTH;


    for n = 1:nActiveReplays % each active replay
        ts = TracesOut.Timestamps{1}(TrialIdxOL(n,1:2)); % start and stop timestamps of the active replay
        % find the relevant indices in the predictions
        [~,idx1] = min(abs(ClosedloopFitTimestamps-ts(1)));
        [~,idx2] = min(abs(ClosedloopFitTimestamps-ts(2)));

        myPSTH   = ClosedloopPSTHs(whichUnit,idx1:idx2); % observed
        myPSTH   = [2 numel(myPSTH) myPSTH];
        ObservedFRs{whichUnit}(n+1,1:numel(myPSTH)) = myPSTH;

        myPSTH   = PredictedClosedloopPSTHs{4}(whichUnit,idx1:idx2); % predicted
        myPSTH   = [2 numel(myPSTH) myPSTH];
        PredictedFRs{whichUnit}(n+1,1:numel(myPSTH)) = myPSTH;

    end

    for n = 1:nReplays % each passive replay
        ts = PassiveOut.Timestamps{1}(ReplayIdx(n,1:2)); % start and stop timestamps of the passive replay
        % find the relevant indices in the prediction
        [~,idx1] = min(abs(PassiveFitTimestamps-ts(1)));
        [~,idx2] = min(abs(PassiveFitTimestamps-ts(2)));

        myPSTH   = PassivePSTHs(whichUnit,idx1:idx2); % observed
        myPSTH   = [3 numel(myPSTH) myPSTH];
        ObservedFRs{whichUnit}(n+1+nActiveReplays,1:numel(myPSTH)) = myPSTH;

        myPSTH   = PredictedPassivePSTHs{4}(whichUnit,idx1:idx2); % predicted
        myPSTH   = [3 numel(myPSTH) myPSTH];
        PredictedFRs{whichUnit}(n+1+nActiveReplays,1:numel(myPSTH)) = myPSTH;
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
multfact = 750/PSTHbinsize;
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
