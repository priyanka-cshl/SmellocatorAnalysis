%% script to analyze basic active/passive replay sessions
SessionName = 'O3_20211005_r0';

%% Housekeeping
Paths = WhichComputer();


%% load the wdw processed sessions

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
    % find the stretches of replay
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

%% plot an example unit
unitIDs = cell2mat({SingleUnits.id})';
%whichUnit = 1;
whichUnit = find(unitIDs==1095);
whichUnit = find(unitIDs==921);
whichUnit = find(unitIDs==1218);
whichUnit = find(unitIDs==979);
thisUnitSpikes = SingleUnits(whichUnit).spikes;
nTrials = nActiveReplays + nReplays + 1;

figure; 
subplot(3,1,1); % plot the observed spikes
% plot the trial structure
PlotBehavior(ReplayTS,[],[],[],[],ReplayTrialVec,[],nTrials);
set(gca,'YLim',[-0.4 nTrials+0.4],'YTick',[],'TickDir','out','XLim',[ReplayTS(1) ReplayTS(end)],'XTick', []);
axis manual;
hold on;

% plot the active replay spikes
for n = 1:nActiveReplays % each passive replay
    ts = TracesOut.Timestamps{1}(TrialIdxOL(n,1:2));
    thisReplaySpikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<=ts(2)));
    thisReplaySpikes = thisUnitSpikes(thisReplaySpikes) - ts(1);
    PlotRaster_v2(thisReplaySpikes,n+1,Plot_Colors('r'));
    %plot(thisReplaySpikes,(n+1)*ones(numel(thisReplaySpikes),1),'*r');
end

% plot the passive replay spikes
for n = 1:nReplays % each passive replay
    ts = PassiveOut.Timestamps{1}(ReplayIdx(n,1:2));
    thisReplaySpikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<=ts(2)));
    thisReplaySpikes = thisUnitSpikes(thisReplaySpikes) - ts(1);
    PlotRaster_v2(thisReplaySpikes,nActiveReplays+n+1,Plot_Colors('t'));
    %plot(thisReplaySpikes,(n+1)*ones(numel(thisReplaySpikes),1),'*r');
end

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
        thisSubtrialSpikes = thisSubtrialSpikes - extraindices(n)*binsize;
    end
    TemplateSpikesAdjusted = vertcat(TemplateSpikesAdjusted,thisSubtrialSpikes);
end
PlotRaster_v2(TemplateSpikesAdjusted,1,Plot_Colors('k'));
%plot(TemplateSpikesAdjusted, ones(numel(TemplateSpikesAdjusted),1),'*k');


% plot the predicted passive replay and closed loop spikes


for whichplot = 2:3
    subplot(3,1,whichplot); % plot the expected spikes
    % plot the trial structure
    PlotBehavior(ReplayTS,[],[],[],[],ReplayTrialVec,[],nTrials);
    set(gca,'YLim',[-0.4 nTrials+0.4],'YTick',[],'TickDir','out','XLim',[ReplayTS(1) ReplayTS(end)],'XTick', []);
    axis manual;
    hold on;
end

% for the active replays - get spike times given the predicted FR from the GLM
for n = 1:nActiveReplays
    ts = TracesOut.Timestamps{1}(TrialIdxOL(n,1:2));
    % find the relevant indices in the predicts
    [~,idx1] = min(abs(ClosedloopFitTimestamps-ts(1)));
    [~,idx2] = min(abs(ClosedloopFitTimestamps-ts(2)));
    %myPSTH   = PassivePSTHs(whichUnit,idx1:idx2);

    subplot(3,1,2);
    myPSTH   = PredictedClosedloopPSTHs{4}(whichUnit,idx1:idx2);
    predictedSpikes = PechePourPoisson(multfact*myPSTH,dt);
    PlotRaster_v2(predictedSpikes,n+1,Plot_Colors('r'));

    subplot(3,1,3);
    myPSTH   = WulfClosedloopPSTHs(whichUnit,idx1:idx2);
    predictedSpikes = PechePourPoisson(multfact*myPSTH,dt);
    PlotRaster_v2(predictedSpikes,n+1,Plot_Colors('r'));
    %plot(predictedSpikes,(n+1+nReplays)*ones(numel(predictedSpikes),1),'*b');
end

% for the passive replays - get spike times given the predicted FR from the GLM
for n = 1:nReplays
    ts = PassiveOut.Timestamps{1}(ReplayIdx(n,1:2));
    % find the relevant indices in the predicts
    [~,idx1] = min(abs(PassiveFitTimestamps-ts(1)));
    [~,idx2] = min(abs(PassiveFitTimestamps-ts(2)));
    %myPSTH   = PassivePSTHs(whichUnit,idx1:idx2);

    subplot(3,1,2);
    myPSTH   = PredictedPassivePSTHs{4}(whichUnit,idx1:idx2);
    predictedSpikes = PechePourPoisson(multfact*myPSTH,dt);
    PlotRaster_v2(predictedSpikes,nActiveReplays+n+1,Plot_Colors('t'));

    subplot(3,1,3);
    myPSTH   = WulfPassivePSTHs(whichUnit,idx1:idx2);
    predictedSpikes = PechePourPoisson(multfact*myPSTH,dt);
    PlotRaster_v2(predictedSpikes,nActiveReplays+n+1,Plot_Colors('t'));
    %plot(predictedSpikes,(n+1+nReplays)*ones(numel(predictedSpikes),1),'*b');
end

for whichplot = 2:3
    % plot the closed loop
    % first predict spike times for the undistorted template stretch
    ts = TracesOut.Timestamps{1}([TemplateIdx(1,3) TemplateIdx(end,2)]); % entire set of template trials
    [~,idx1] = min(abs(ClosedloopFitTimestamps-ts(1)));
    [~,idx2] = min(abs(ClosedloopFitTimestamps-ts(2)));

    subplot(3,1,whichplot);
    switch whichplot
        case 2
            myPSTH = PredictedClosedloopPSTHs{4}(whichUnit,idx1:idx2);
        case 3
            myPSTH = WulfClosedloopPSTHs(whichUnit,idx1:idx2);
    end
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
            thisSubtrialSpikes = thisSubtrialSpikes - extraindices(n)*binsize;
        end
        PredictedTemplateSpikesAdjusted = vertcat(PredictedTemplateSpikesAdjusted,thisSubtrialSpikes);
    end
    PlotRaster_v2(PredictedTemplateSpikesAdjusted,1,Plot_Colors('k'));
    %plot(TemplateSpikesAdjusted, ones(numel(TemplateSpikesAdjusted),1),'*k');
end

