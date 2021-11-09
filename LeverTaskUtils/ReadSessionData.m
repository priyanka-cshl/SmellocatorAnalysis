% test script to extract behavior data and replot session
% for animals trained on the fixed gain version of the task (post 08.2018)

function [MyData, MyParams, DataTags] = ReadSessionData(FileName, PIDflag)
if nargin<2
    PIDflag = 0;
end

[MyData, MyParams, DataTags, TFtype] = LoadSessionData(FileName, 0, PIDflag);

%% clean up params table
% HACK 1: only keep entries that have non-zero timestamps
MyParams(1: find(MyParams(:,1)==0,1,'last')-1,:) = [];
% HACK 2: params are written to with -ve timestamp before the 
% actual update happens from the Arduino
% if timestamps are negative - ignore those
MyParams(find(MyParams(:,1)<0),:) = [];

%% HACK: for first day of block shifts - fix target zone definitions
f = find(MyParams(:,23)~=121);
if ~isempty(f)
    U = unique(MyParams(:,18:20),'rows');
    for i = 1:size(U,1)
        if rem(U(i,2),0.25)
            x = find(MyParams(:,19) == U(i,2));
            y = find(abs(MyParams(:,19) - U(i,2) + 0.6)<=0.1); 
            for j = 1:numel(x)
                MyParams(x(j),18:20) = MyParams(y(1),18:20); 
            end
        end
    end
end

%% open loop and replays
OL_Blocks  = numel(find(diff(MyParams(:,32))==1)); % in timestamps
OL_template = []; 
Replay_Trials = [];
if any(OL_Blocks)
    for OL = 1:OL_Blocks
        foo = [];
        foo(:,1) = (find(diff(MyParams(:,32))== 1)+1):(find(diff(MyParams(:,32))== -1));
        foo(:,2) = OL;
        OL_template = vertcat(foo,OL_template);
    end
    
    % any replays?
    Replay_Trials = find(diff(MyParams(:,32))== 2) + 1;
end

%% for each trial
for thisTrial = 1:size(MyParams,1)
    
    %% Fill targetzone values into col2,3 - for plotting behavior GUI style
    if thisTrial < size(MyParams,1)
        f = find((MyData(:,1)>=MyParams(thisTrial,1)) &...
            (MyData(:,1)<MyParams(thisTrial+1,1)));
    else
        f = find(MyData(:,1)>=MyParams(end,1));
    end
    MyData(f,2) = MyParams(thisTrial,18);
    MyData(f,3) = MyParams(thisTrial,20);
    
    %% detect and categorize perturbations
    if MyParams(thisTrial,26)>1 % was a perturbed trial
        switch MyParams(thisTrial,26)
            case 2 % fake zone
                MyData(f,11) = MyParams(thisTrial,28); %FZoneHighLim
                MyData(f,12) = MyParams(thisTrial,30); %FZoneLowLim
            case 3 % to detect NoOdor trials
                MyData(f,11) = 100*MyParams(thisTrial,26); %300
                MyData(f,12) = -1;
            case 4 % flip mapping
                MyData(f,11) = 100*MyParams(thisTrial,26); %400
                MyData(f,12) = -1;
            case {5,6,7} % location offset
                MyData(f,11) = 100*MyParams(thisTrial,26);
                MyData(f,12) = MyParams(thisTrial,27); % offset size
            case 8 % gain change
                MyData(f,11) = 100*MyParams(thisTrial,26);
                MyData(f,12) = MyParams(thisTrial,27); % gain value
            case 9 % feedback halt
                %MyData(f,11) = 100*MyParams(thisTrial,26);
                MyData(f,11) = 6; %FZoneHighLim
                MyData(f,12) = 5; %FZoneLowLim
            case 10 % feedback pause or feedback halt (2021 or later)
                %MyData(f,11) = 100*MyParams(thisTrial,26);
                if str2double(FileName(regexp(FileName,'_','once')+(1:4))) > 2020
                    % feedback halt perturbation with location flip
                    MyData(f,11) = 100*MyParams(thisTrial,26); %1000
                    MyData(f,12) = MyParams(thisTrial,28) - 121; % halted odor location (w.r.t. center)
                else
                    % feedback pause perturbation
                    MyData(f,11) = 6; %FZoneHighLim
                    MyData(f,12) = 6; %FZoneLowLim
                end
            case 13 % LED only trials
                MyData(f,11) = 100*MyParams(thisTrial,26); % 1300
                
        end
    end
    
    % flag block shift trials
    if MyParams(thisTrial,23)~=121 % was a perturbed trial
        MyData(f,11) = 100*11;
        MyData(f,12) = MyParams(thisTrial,23)-121; % shift size
    end
    
    % detect rule reversal trials
    TrialState = MyData(f,find(ismember(DataTags,'TrialON')));
    MotorState = MyData(f,find(ismember(DataTags,'Motor')));
    OdorStartState = MotorState(find(diff(TrialState>0),1,'first')+1);
    if ~isempty(OdorStartState)
        OdorStartState = (OdorStartState>0); % this is true if TFtype = 1
        if OdorStartState ~= TFtype
            MyData(f,11) = 1400;
            MyData(f,12) = OdorStartState+1;
        end
    end
    
    % flag open loop template trials
    if ismember(thisTrial,OL_template(:,1))
        MyData(f,11) = 1500;
        MyData(f,12) = -1;
    end
    
    % flag replay trials
    if ismember(thisTrial,Replay_Trials)
        MyData(f,11) = 1600;
        MyData(f,12) = -2;
    end
    
end

%% convert trial_ON column to odor IDs
% column number = 6 in MyData
for odor = 1:4
    f = find((MyData(:,6)>=odor^2) & (MyData(:,6)<(odor+1)^2));
    MyData(f,6) = odor;
end

%% HACK: to compensate for code bug for No odor trials - old - not currently used I think - PG (2021-10-07)
if isempty(strfind(FileName,'LR')) % skip this step for visual version
    if any(MyParams(:,2)==0) && ~any(MyData(:,6)>=4)
        % cheat to prevent last trial from being NoOdor Trial
        if MyParams(end,2) == 0
            MyParams(end,:) = [];
        end
        % extra hack to figure out when the odor was on if odor ID = 0
        NoOdorTrials(:,1) = MyParams(find(MyParams(:,2)==0),1);
        NoOdorTrials(:,2) = MyParams(find(MyParams(:,2)==0)+1,1);
        NoOdorTrials(:,3) = find(MyParams(:,2)==0);
        
        if size(NoOdorTrials,1)>0
            for thisTrial = 1:size(NoOdorTrials,1)
                indices = find((MyData(:,1)>NoOdorTrials(thisTrial,1)) & (MyData(:,1)<NoOdorTrials(thisTrial,2)));
                trialthreshold = MyParams(NoOdorTrials(thisTrial,3),11);
                trialhold = round(MyParams(NoOdorTrials(thisTrial,3),13))/2; % convert to samples
                thistriallever = MyData(indices,4);
                thistriallever(thistriallever<trialthreshold) = 0;
                thistriallever(thistriallever>=trialthreshold) = 1;
                triggerstart = find(diff(thistriallever)==1);
                triggerstop = find(diff(thistriallever)==-1);
                m = 1;
                while m<=min([numel(triggerstart),numel(triggerstop)])
                    if (triggerstop(m)-triggerstart(m)+1)>=trialhold
                        indices(1:triggerstop(m)-1,:) = [];
                        MyData(indices,6) = 4; % odor 4
                        break;
                    end
                    m = m + 1;
                end
            end
        end
    end
end

end