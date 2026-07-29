%% script to plot PID traces for CID sessions done on photoncerber 

%recDir = '/mnt/grid-hs/mdussauz/CID/PID/2022-05-31_08-54-29_16Odor/Record Node 104';
recDir = '/mnt/grid-hs/mdussauz/CID/PID/newCID/2025-05-14_12-10-02/Record Node 104';
%recDir = '/mnt/grid-hs/mdussauz/CID/PID/newCID/2025-05-14_15-20-57/Record Node 104';

% load the preprocessed data (or preprocess first)
[PIDTrace, TTLs, StimSettings] = getPID_CID(recDir);

%% organize by trial, stimulus etc

% some settings
SamplingRate = 500; % in Hz

switch StimSettings.SessionType
    case '16_Odors'
        nStim = size(StimSettings.Odors,1); % no. of odors delivered
        OdorList = StimSettings.Odors;
        nTypes = numel(StimSettings.Dilutions); % concentrations used
        ConcList = StimSettings.Dilutions;
        window = [-StimSettings.timing(2) sum(StimSettings.timing(4:5))]/1000; % in seconds [-preOdor odor+post]
    case 'newCID'
        nStim = size(StimSettings.miniOdors,1) + size(StimSettings.megaOdors,1); % no. of odors delivered
        OdorList = [StimSettings.miniOdors; StimSettings.megaOdors];
        nTypes = numel(StimSettings.Dilutions); % concentrations used
        ConcList = StimSettings.Dilutions;
        window = [-StimSettings.timing(2) sum(StimSettings.timing(4:5))]/1000; % in seconds [-preOdor odor+post]
    case 'ConcentrationSeries'
        window = [-StimSettings.timing(2) sum(StimSettings.timing(4:5))]/1000; % in seconds [-preOdor odor+post]
end

window = window -1; % 1 sec buffer
baselineWindow = 1*SamplingRate; % samples
PIDOut = []; TimeTraceOut = []; SteadyState = []; MeanVals = []; TimeOut = []; prevOdor = [];

for odor = 1:nStim % every odor
    for conc = 1:nTypes % at every conc.
        whichTrials = find( (TTLs.Trial(:,4) == OdorList(odor)) & (TTLs.Trial(:,5) == ConcList(conc)) );
        for n = 1:numel(whichTrials)
            t = TTLs.Trial(whichTrials(n),7:8); % odor start, odor stop
            t = t + window;
            [~,idx1] = min(abs(PIDTrace(:,1)-t(1)));
            [~,idx2] = min(abs(PIDTrace(:,1)-t(2)));
            thisTrialPID = PIDTrace(idx1:idx2,2);

            PIDOut(n,1:numel(thisTrialPID),odor,conc) = thisTrialPID;
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

            % preodor baseline
            [~,idx3] = min(abs(PIDTrace(:,1)-TTLs.Trial(whichTrials(n),7)));
            baselines(n,1,odor,conc) = median(PIDTrace(idx3+[-baselineWindow 0],2));
        end
    end
end

%% Plotting related settings

if strcmp(StimSettings.SessionType,'ConcentrationSeries')
    stackConcs = 1;
    if stackConcs == 2
        newStimCol = (TTLs.Trial(:,5)*10^4) + TTLs.Trial(:,4);
        newTypesCol = TTLs.Trial(:,5)*0;
        TTLs.Trial(:,4:5) = [newStimCol newTypesCol];
        StimSettings.SessionType = '16_Concs*';
    end
    if stackConcs == 1
        newStimCol = TTLs.Trial(:,5) + TTLs.Trial(:,4);
        newTypesCol = TTLs.Trial(:,5)*0;
        TTLs.Trial(:,4:5) = [newStimCol newTypesCol];
        StimSettings.SessionType = '16_Concs';
    end
end

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

%% plotting
trialgradient = 0;
normalized = 0;
baselinesubtraced = 1;
if ~baselinesubtraced
    normalized = 0;
end
AverageOut = [];
traceLength = (diff(window)+(StimSettings.timing(3)/1000))*SamplingRate;
% Plot the raw values, later trials are darker
figure;
for odor = 1:nStim
    subplot(4,4,odor);
    hold on
    YVals = [];
    for conc = 1:nTypes
        AverageOut(odor,:,conc) = zeros(traceLength,1);
        reps = 0;
        for n = 1:size(PIDOut,1)
            if baselines(n,1,odor,conc)
                thisTrialPID = [];
                thisTrialPID(:,2) = squeeze(PIDOut(n,:,odor,conc));
                thisTrialPID(1:numel(TimeOut),1) = TimeOut;
                untilIdx = find(squeeze(TimeTraceOut(n,:,odor,conc))==0,1,'first');
                thisTrialPID(untilIdx:end,:) = [];
                if baselinesubtraced
                    thisTrialPID(:,2) = thisTrialPID(:,2) - baselines(n,1,odor,conc); % subtract baseline
                end
                if normalized
                    thisTrialPID(:,2) = thisTrialPID(:,2)/max(thisTrialPID(:,2));
                end
                thisTrialPID(1:numel(TimeOut),1) = TimeOut;
                thisTrialPID(numel(TimeOut)+1:end,:) = [];
                if trialgradient
                    thisTrialcolor = mycolors(odor,:) + (1 - mycolors(odor,:))*(1-n*0.2);
                    plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',thisTrialcolor);
                else
                    plot(thisTrialPID(:,1),thisTrialPID(:,2),'Color',mycolors(odor,:));
                end
                YVals(n,:) = [min(thisTrialPID(:,2)) max(thisTrialPID(:,2))];
                AverageOut(odor,:,conc) = AverageOut(odor,:,conc) + thisTrialPID(1:traceLength,2)';
                reps = reps+1;
            end
        end
        AverageOut(odor,:,conc) = AverageOut(odor,:,conc)/reps;
        %plot(thisTrialPID(1:traceLength,1),squeeze(AverageOut(odor,:,conc)),'Color',mycolors(odor,:),'LineWidth',2);
        plot(thisTrialPID(1:traceLength,1),squeeze(AverageOut(odor,:,conc)),'Color','k');
    end
    YVals(n+1,:) = [min(YVals(:,1)) max(YVals(:,2))];
    YVals(n+1,:) = YVals(n+1,:) + (abs(diff(YVals(n+1,:)))*0.1)*[-1 1];
    set(gca,'XLim',[-10 30],'YLim',YVals(n+1,:));
end

%%
skip = 1
if ~skip
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
end