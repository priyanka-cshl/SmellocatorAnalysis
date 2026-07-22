%% script to plot PID traces for CID sessions done on photoncerber 

recDir = '/mnt/grid-hs/mdussauz/PID_test/2022-05-31_08-54-29_16Odor/Record Node 104';
%recDir = '/mnt/grid-hs/mdussauz/PID_test/2022-05-30_18-19-56/Record Node 104'
[PIDTrace, TTLs, StimSettings] = getPID_CID(recDir);

if strcmp(StimSettings.SessionType,'ConcentrationSeries')
stackConcs = 1;
if stackConcs == 2
    if strcmp(StimSettings.SessionType,'ConcentrationSeries')
        newStimCol = (TTLs.Trial(:,5)*10^4) + TTLs.Trial(:,4);
        newTypesCol = TTLs.Trial(:,5)*0;
        TTLs.Trial(:,4:5) = [newStimCol newTypesCol];
        StimSettings.SessionType = '16_Concs*';
    end
end
if stackConcs == 1
    if strcmp(StimSettings.SessionType,'ConcentrationSeries')
        newStimCol = TTLs.Trial(:,5) + TTLs.Trial(:,4);
        newTypesCol = TTLs.Trial(:,5)*0;
        TTLs.Trial(:,4:5) = [newStimCol newTypesCol];
        StimSettings.SessionType = '16_Concs';
    end
end
end

nStim = numel(unique(TTLs.Trial(:,4))); % no. of odors delivered
OdorList = unique(TTLs.Trial(:,4));
nTypes = numel(unique(TTLs.Trial(:,5))); % concentrations used
ConcList = unique(TTLs.Trial(:,5));

switch StimSettings.SessionType
    case {'newCID', '16_Odors'}
        mycolors = brewermap(nStim,'Dark2');
        foo = round(nStim/2)+1;
        mycolors(foo:end,:) = mycolors(foo:end,:)/2;

    case {'16_Concs*'}
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        for c = 0:0.2:0.6
            lighter = basecolors + (1 - basecolors)*c;
            mycolors = vertcat(lighter, mycolors);
        end
        StimSettings.SessionType = '16_Odors';

    case '16_Concs'
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        scale = 0.6:-0.2:0;
        for c = 1:4
            mycolors(:,:,c) = basecolors + (1 - basecolors)*scale(c);
        end
        mycolors = reshape(permute(mycolors, [3 1 2]), 20, 3);      % rows: [color1's N shades; color2's N shades; ...]
        StimSettings.SessionType = '16_Odors';
end

%%
window = [-StimSettings.timing(2) sum(StimSettings.timing(4:5))]/1000; % in seconds
window = window -1;
PIDOut = []; TimeTraceOut = []; SteadyState = []; MeanVals = []; TimeOut = []; prevOdor = [];
for odor = 1:nStim
    for conc = 1:nTypes
        whichTrials    = find( (TTLs.Trial(:,4) == OdorList(odor)) & (TTLs.Trial(:,5) == ConcList(conc)) );
        for n = 1:numel(whichTrials)
            t = TTLs.Trial(whichTrials(n),7:8); % odor start, odor stop
            t = t + window;
            [~,idx1] = min(abs(PIDTrace(:,1)-t(1)));
            [~,idx2] = min(abs(PIDTrace(:,1)-t(2)));
            thisTrialPID = PIDTrace(idx1:idx2,2);
            PIDOut(n,1:numel(thisTrialPID),odor,conc) = thisTrialPID;
%             PIDvals{n,odor+1}  = PIDTrace(idx1:idx2,2);
%             SteadyState(n,odor+1) = mean(PIDvals{n,odor+1}(end-99:end));
            TimeTraceOut(n,1:numel(thisTrialPID),odor,conc) = PIDTrace(idx1:idx2,1);
            if (whichTrials(n)-1)>0
                prevOdor(n,1,odor,conc) = TTLs.Trial(whichTrials(n)-1,4);
            else
                prevOdor(n,1,odor,conc) = TTLs.Trial(whichTrials(n),4);
            end
            if odor == 1 && conc == 1 && n == 1
                TimeOut = PIDTrace(idx1:idx2,1);
                TimeOut = TimeOut - TTLs.Trial(whichTrials(n),7);
            end
        end
    end
end

%% plotting
%%
% figure;
% conc = 1;
% for odor = 1:nStim
%     subplot(4,4,odor);
%     hold on
%     for n = 1:size(PIDOut,1)
%         thisTrialPID = [];
%         thisTrialPID(:,2) = squeeze(PIDOut(n,:,odor,conc));
%         thisTrialPID(1:numel(TimeOut),1) = TimeOut;
%         thisTrialPID(numel(TimeOut)+1:end,:) = []; 
%         plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',mycolors(odor,:));
%     end
% end

figure;
conc = 1;
for odor = 1:nStim
    subplot(4,4,odor);
    hold on
    YVals = [];
    for n = 1:size(PIDOut,1)
        thisTrialPID = [];
        thisTrialPID(:,2) = squeeze(PIDOut(n,:,odor,conc));
        thisTrialPID(1:numel(TimeOut),1) = TimeOut;
        thisTrialPID(numel(TimeOut)+1:end,:) = []; 
        thisTrialcolor = mycolors(odor,:) + (1 - mycolors(odor,:))*(1-n*0.2);
        %plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',mycolors(odor,:));
        plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',thisTrialcolor);
        YVals(n,:) = [min(thisTrialPID(:,2)) max(thisTrialPID(:,2))];
    end
    YVals(n+1,:) = [min(YVals(:,1)) max(YVals(:,2))];
    YVals(n+1,:) = YVals(n+1,:) + (abs(diff(YVals(n+1,:)))*0.1)*[-1 1];
    set(gca,'XLim',[-10 30],'YLim',YVals(n+1,:));
end

deltaT = TimeOut(end) - TimeOut(1) + 0.1;
figure;
conc = 1;
for odor = 1:nStim
    subplot(4,4,odor);
    hold on
    for n = 1:size(PIDOut,1)
        thisTrialPID = [];
        thisTrialPID(:,2) = squeeze(PIDOut(n,:,odor,conc));
        thisTrialPID(1:numel(TimeOut),1) = TimeOut + deltaT*(n-1);
        thisTrialPID(numel(TimeOut)+1:end,:) = []; 
        thisTrialcolor = mycolors(odor,:) + (1 - mycolors(odor,:))*(1-n*0.2);
        plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',mycolors(odor,:));
        %plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',thisTrialcolor);
    end
end

figure;
conc = 1;
for odor = 1:nStim
    subplot(4,4,odor);
    hold on
    for n = 1:size(PIDOut,1)
        thisTrialPID = [];
        thisTrialPID(:,2) = squeeze(PIDOut(n,:,odor,conc));
        thisTrialPID(1:numel(TimeOut),1) = TimeOut + deltaT*(n-1);
        thisTrialPID(numel(TimeOut)+1:end,:) = []; 
        thisTrialcolor = mycolors(squeeze(prevOdor(n,1,odor,conc)),:);
        %plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',mycolors(odor,:));
        plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',thisTrialcolor);
    end
end

figure;
conc = 1;
for odor = 1:nStim
    %subplot(4,4,odor);
    hold on
    for n = 1:size(PIDOut,1)
        thisTrialPID = [];
        thisTrialPID(:,2) = squeeze(PIDOut(n,:,odor,conc));
        thisTrialPID(:,1) = squeeze(TimeTraceOut(n,:,odor,conc));
        thisTrialPID(thisTrialPID(:,1)==0,:) = [];
        %thisTrialPID(numel(TimeOut)+1:end,:) = []; 
        thisTrialcolor = mycolors(odor,:) + (1 - mycolors(odor,:))*(1-n*0.2);
        plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',mycolors(odor,:));
        %plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',thisTrialcolor);
    end
end

%%
figure;
conc = 1;
for odor = 1:nStim
    subplot(1,5,floor(OdorList(odor)));
    hold on
    for n = 1:size(PIDOut,1)
        thisTrialPID = [];
        thisTrialPID(:,2) = squeeze(PIDOut(n,:,odor,conc));
        thisTrialPID(1:numel(TimeOut),1) = TimeOut;
        thisTrialPID(numel(TimeOut)+1:end,:) = []; 
        plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',mycolors(odor,:));
    end
end

% count = 0;
% for loc = -90:15:90
%     count = count + 1;
%     subplot(1,13,count);
%     hold on;
%     f = find(whichlocations==loc);
%     for i = 1:numel(f)
%         plot(PIDvals{f(i),odor+1});
%     end
%     set(gca,'YLim',[-1.5 1]);
%     MeanVals{odor+1}(count,1) = mean(SteadyState(f,odor+1));
%     MeanVals{odor+1}(count,2) = std(SteadyState(f,odor+1));
% end
    