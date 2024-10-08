function [OpenLoop] = ExtractReplayTrials_v2(Traces, TrialInfo, TTLs, ReplayTTLs, PassiveTraces)

%% globals
global SampleRate; % = 500; % samples/second
global startoffset; % = 1; % seconds
traceOverlap = SampleRate*startoffset;
AllTargets = 1:0.25:3.75; % for assigning target zone values

% rename all trials flagged as Halt-Flip-Templates or Offset-II-templates as OL-Template
PerturbedTemps = [];
if any(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
    PerturbedTemps = find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip-Template'));
elseif any(strcmp(TrialInfo.Perturbation(:,1),'Offset-II-Template'))
    PerturbedTemps = find(strcmp(TrialInfo.Perturbation(:,1),'Offset-II-Template'));
end
for i = 1:numel(PerturbedTemps)
    TrialInfo.Perturbation{PerturbedTemps(i),1} = 'OL-Template';
end

% sometimes odorstart is computed as NaN - mouse started moving too quickly
% force it to be zero
TrialInfo.OdorStart(find(isnan(TrialInfo.OdorStart(:,1))),1) = 0;

% How many open loop templates are there?
% typically only 1 - in rare cases there might be two
TemplateTrials(:,1) = find(diff(strcmp(TrialInfo.Perturbation,'OL-Template'))== 1) + 1;
if ~strcmp(TrialInfo.Perturbation(end,1),'OL-Template')
    TemplateTrials(:,2) = find(diff(strcmp(TrialInfo.Perturbation,'OL-Template'))==-1);
else
    TemplateTrials(:,2) = [find(diff(strcmp(TrialInfo.Perturbation,'OL-Template'))==-1); size(TrialInfo.Perturbation,1)];
end

% for each template - get a concatenated trace for
% Lever, Motor, Resp, Licks, TrialON, Rewards, TargetZone, Timestamps (?)
whichTraces = fieldnames(Traces);
%for i = 1:size(whichTrials) % no. of passive replays
%% passive replays
if ~exist('TrialOFF','var') % no active replays
    if ~isempty(PassiveTraces)
        OpenLoop.PassiveReplayTraces = PassiveTraces;
    end

    whichTrials = ReplayTTLs.TrialID;
    whichTrials(:,whichTrials<size(TrialInfo.Odor,1)) = [];

    TrialOFF = [];
    Nsubtrials = mode(diff(TemplateTrials',1)) + 1;
    for r = 1:numel(whichTrials) % for every replay
        OpenLoop.PassiveReplayTraces.TrialIDs{r} = whichTrials(r);
        % the replay trace is already concatenated - but need to know
        % how to split it into subtrials such that the template can
        % be aligned to it
        % this is because a few samples get missed during Arduino
        % updates at the end of every OL - template trial;

        % 1. get Trial OFF indices w.r.t. to trace start (includes startoffset)
        if isempty(ReplayTTLs) % Use the behavior trace itself
            % find reward timestamps for the replayed trace
            % use reward TS as proxy of trial OFF
            TrialOFF(:,r) = find(diff(Traces.Rewards{whichTrials(r)})==1)+1;
        else % use the OdorOFF TTLs from OEPS as proxy for trial OFF
            % which set of ReplayTTLs to use
            f = find(ReplayTTLs.TrialID==whichTrials(r));
            T_off = ReplayTTLs.OdorValve{f}(:,2);
            % sometimes there is 1 extra odor transition at trialstart
            % if replay trial starts with a different odor - delete that
            while size(T_off,1)>Nsubtrials
                T_off(1,:) = [];
            end
            % account for the padded startoffset samples in the template
            TrialOFF(:,r) = ceil((T_off+startoffset)*SampleRate);

            OpenLoop.PassiveReplayTraces.OdorSeq(r,1:size(ReplayTTLs.OdorValve{f},1)) = ReplayTTLs.OdorValve{f}(:,4);
        end
    end
end
%end
   
for i = 1:size(TemplateTrials,1) % no. of templates
    %% Open Loop Template
    % get trial IDs for trials that constitute the replay stretch
    whichTrials = TemplateTrials(i,1):TemplateTrials(i,2);
    OpenLoop.TemplateTraces.TrialIDs{i} = whichTrials;
    
    % get all traces and concatenate them
    for j = 1:size(whichTraces,1)
        temp = cellfun(@(x) ...
            x(1:end-traceOverlap), Traces.(whichTraces{j})(whichTrials), ...
            'UniformOutput', false);
        
        % add in the overlap for the very last trial if needed
        % also make sure that there are atleast startoffset*SampleRate
        % samples after TrialOFF
        OpenLoop.TemplateTraces.(whichTraces{j})(i) = {[cell2mat(temp(:)); ...
            Traces.(whichTraces{j}){whichTrials(end)}(end-traceOverlap+1:end,1)]};
        
    end
    % use the TrialON column
    % 1. to construct the targetzone trace - for plotting
    % 2. to include OdorON periods (in -ve)
    TrialTrace = cell2mat(OpenLoop.TemplateTraces.Trial(i));
    % make sure any trial ON periods preceding trial1 start are ignored
    TrialTrace(1:traceOverlap,1) = 0;
    % get trial ON-OFF indices
    Idx = [find(diff(TrialTrace>0)==1)+1 find(diff(TrialTrace>0)==-1)];
    % get OdorStart Times w.r.t. Trial start (from the behavior file) 
    Idx(:,3) = Idx(:,1) + ceil(SampleRate*TrialInfo.OdorStart(whichTrials,1));
    if ~isempty(TTLs)
        % OdorStart Times can also be extracted from OEPS TTLs
        Idx(:,4) = Idx(:,1) + ceil(SampleRate*TTLs.Trial(whichTrials,4));
        
        % If Timestamps were dropped - replace odor start values with those
        % from the Ephys file
        if any(TrialInfo.TimeStampsDropped(whichTrials))
            Idx(find(TrialInfo.TimeStampsDropped(whichTrials)),3) = Idx(find(TrialInfo.TimeStampsDropped(whichTrials)),4);
        end
        if numel(find(abs(Idx(:,3)-Idx(:,4))>5))>1
            disp('mismatch in OdorOn Timestamps between behavior and OEPS files');
            %keyboard;
        end
    end
    % sometimes if the mouse takes too long to start the first trial (>1sec), so
    % odorON indices can be negative
    Idx(1,find(Idx(1,:)<0)) = 1;
    temp = 0*TrialTrace;
    x1 = 1;
    
    for k = 1:size(Idx,1) % every subtrial
        % TargetZone vector
        x2 = Idx(k,2);
        temp(x1:x2,1) = AllTargets(TrialInfo.TargetZoneType(whichTrials(k)));
        x1 = x2 + 1;
        
        % OdorON + TrialON vector
        TrialTrace(Idx(k,3):Idx(k,1)-1,1) = -TrialTrace(Idx(k,1));
    end
    OpenLoop.TemplateTraces.TargetZone(i) = {temp};
    OpenLoop.TemplateTraces.Trial(i) = {TrialTrace};
    
    %% Aligning open loop template and replay traces
    % if there are any corrupt replays we want to ignore them when deciding
    % where to chop the template for extra samples
    X = [];
    if isempty(ReplayTTLs)
        X(:,2) = find(diff(OpenLoop.TemplateTraces.Rewards{i})==1)+1;
    else
        X(:,2) = Idx(:,2);
    end
    whichOdors = TrialInfo.Odor(whichTrials);
    whichPassive = find(ismember(OpenLoop.PassiveReplayTraces.OdorSeq(:,1:numel(whichOdors)),whichOdors','rows'));
    if exist('TrialOFF','var')
        %X(:,4) = mode(TrialOFF(1,:)) + mode(TrialOFF - TrialOFF(1,:),2);
        if ~isempty(whichPassive)
            X(:,4) = TrialOFF(:,whichPassive); % + TrialOFF(:,whichPassive) - TrialOFF(1,:);
        end
    else
        X(:,4) = X(:,2); % for sessions where no active replays were done
    end
    % col 2 and col 4 now contain stop indices for each subtrial
    % col 2 - w.r.t. template, col 4 - w.r.t. replayed traces
    % get start indices for both sets of subtrials
    X(:,[1 3]) = [[1 1]; 1+X(1:end-1,[2 4])];
    % get subtrial lenths
    X(:,5) = X(:,2) - X(:,1); % template
    X(:,6) = X(:,4) - X(:,3); % replays
%     % sometimes the first replay subtrial has a few extra samples
%     if (X(1,6)-X(1,5))>0
%         X(1,3) = X(1,3) + X(1,6) - X(1,5);
%     end
    % get the count of extra samples in the template - difference in subtrial lengths
    X(:,7) = X(:,5) - X(:,6);
    
    % Replace extra samples by NaNs in the template
    for j = 1:size(whichTraces,1)
        for k = 1:size(X,1) % for every subtrial
            if X(k,7) > 0
                x1 = X(k,1); x2 = X(k,1) - 1 + X(k,7);
                OpenLoop.TemplateTraces.(whichTraces{j}){i}(x1:x2,1) = NaN;
            end
        end
        
        if X(1,7) < 0 % if the first trial has less samples than the replay
            OpenLoop.TemplateTraces.(whichTraces{j}){i} = vertcat( ...
                NaN*ones(X(1,7),1), ...
                OpenLoop.TemplateTraces.(whichTraces{j}){i});
        end
    end
    
    %% sanity checks
    if isfield(OpenLoop,'PassiveReplayTraces')
        foo = OpenLoop.TemplateTraces.Motor{i};
        foo = foo(~isnan(foo));
        % if replay traces are longer than the template traces after accounting
        % for NaNs - delete the extra samples
        
        if length(OpenLoop.PassiveReplayTraces.Motor{whichPassive})>length(foo)
            samps_to_delete = length(OpenLoop.PassiveReplayTraces.Motor{whichPassive})-length(foo) - 1;
            for j = 1:size(whichTraces,1)
                if isfield(OpenLoop.ReplayTraces,whichTraces{j})
                    OpenLoop.ReplayTraces.(whichTraces{j}){i}(end-samps_to_delete:end,:) = [];
                end
            end
        end
        
        x = length(OpenLoop.PassiveReplayTraces.Motor{whichPassive});
        Residuals = OpenLoop.PassiveReplayTraces.Motor{whichPassive} - ...
            foo(1:x,1);
        % ignore residuals before first trial start and last trial end
        Residuals(1:2*(SampleRate*startoffset),:) = 0; % being generous about errors at the beginning
        Residuals(X(end,4):end,:) = 0;

        %ErrorDist = fitdist(Residuals(:),'normal');
        ErrorDist = fitdist(Residuals(~isinf(Residuals(:))),'normal'); % if timestamps were dropped - traces have Infs in them
        % check if mean is ~0 and if sigma is very small (<5)
        if ErrorDist.sigma<5
        %if (ErrorDist.mean<0.1) && ErrorDist.sigma<5
            disp('template and replay traces align well');
        else
            disp('template and replay traces do not seem to align well');
            keyboard;
        end
    end
    
end
OpenLoop.TTLs = ReplayTTLs;
    