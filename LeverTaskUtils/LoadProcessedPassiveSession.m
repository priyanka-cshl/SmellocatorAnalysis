%% Passive replays
if ~isempty(PassiveReplayTraces)
    nReps = size(PassiveReplayTraces.Lever,2);
    
    % Stitch all reps together into one long vector
    whichTraces = fieldnames(PassiveReplayTraces);

    for j = 1:size(whichTraces,1)
        temp = cellfun(@(x) ...
            x(1:end), PassiveReplayTraces.(whichTraces{j})(1:nReps), ...
            'UniformOutput', false);
        
        % add in the overlap for the very last trial if needed
        % also make sure that there are atleast startoffset*SampleRate
        % samples after TrialOFF
        PassiveTracesOut(:,j) = cell2mat(temp(:));
    end
    
    % filter the sniff trace
    PassiveTracesOut(:,3) = filtfilt(b,a,PassiveTracesOut(:,3)); %apply the filter to x(t)
    
    % separate the trial vector into odor periods
    
    % get the template trace
    replayTrials = OpenLoop.TemplateTraces.Trial{1};
    replayTrials(isnan(replayTrials)) = [];
    % only keep trace until last trial was off
    tracelength = find(diff(abs(replayTrials)==-1),1,'last');
    replayTrials(tracelength:end) = [];
    % remove initial samples - before trial ON
    replayTrials(1:find(replayTrials>0,1,'first')-1) = [];
    
    % get start stop times of passive replay in the concat trace
    StartStopIdx = [find(diff(PassiveTracesOut(:,5))<0) find(diff(PassiveTracesOut(:,5))>0)];
    PassiveTrialTrace = PassiveTracesOut(:,5);
    for i = 1:size(StartStopIdx,1)
        x1 = StartStopIdx(i,1);
        x2 = x1 + numel(replayTrials) - 1;
        PassiveTrialTrace(x1:x2) = replayTrials*10;
    end
    PassiveTrialTrace(size(PassiveTracesOut,1)+1:end) = [];
    
    %% Split the trial vector into the three odors
    OdorTrace = PassiveTrialTrace;
    OdorTrace(abs(OdorTrace)>5) = OdorTrace(abs(OdorTrace)>5)/10;
    for i = 1:3
        Motor = PassiveTracesOut(:,find(strcmp(whichTraces,'Motor')))';
        Motor(abs(OdorTrace)~=i) = NaN;
        Motor = 125 - Motor;
        %Motor(isnan(Motor)) = 0;
        PassiveTracesOut(:,end+1) = Motor;
    end
    % during air periods
    Motor = PassiveTracesOut(:,find(strcmp(ColNames,'Motor')))';
    Motor(abs(OdorTrace)~=0) = NaN;
    Motor = 125 - Motor;
    %Motor(isnan(Motor)) = 0;
    PassiveTracesOut(:,end+1) = Motor;

end