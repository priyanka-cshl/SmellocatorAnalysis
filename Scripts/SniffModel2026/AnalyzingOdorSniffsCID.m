% run after running CID Response Prepper 
%[StimSettings, TTLs, SingleUnits, AllSpikes, TrialWiseSniffs, SniffsPlot] = CIDResponsePrepper(myKsDir);

% TTLs.Trial : Columns
% 1 to 2: Trial Start, Stop
% 3     : Trial Duration?
% 4     : Stimulus identity (or mixed if conc. series)
% 5     : Stimulus intensity
% 6     : repeat index
% 7 to 8: Odor On, Off
% 9, 10 : Second Odor Pulse On, Off (new CID)
% 11    : true time of first (or nth) sniff after odor ON
% 12    : relative time of first (or nth) sniff after odor ON

% TrialWiseSniffs: Columns
% 1 to 3: inh start, end, next w.r.t. this Trial's odor Onset
% 4     : trial index
% 5     : sniff index within a trial, 0 = first sniff after odor onset
% 6     : actual odor ON time
% 7     : trial phase : 0 - pre-odor, 1 - odor, 1.5 - second odor pulse, 
%                       2 - post-odor, 2.5 - post second odor pulse

%%
nStim = unique(TTLs.Trial(:,4));
mycolors = brewermap(numel(nStim),'Dark2');
foo = round(numel(nStim)/2)+1;
mycolors(foo:end,:) = mycolors(foo:end,:)/2;
medianSniff = median(TrialWiseSniffs(:,3)-TrialWiseSniffs(:,1));
TrialWiseSniffs(:,9) = TTLs.Trial(TrialWiseSniffs(:,4),4); % stimulus identity
odorList = [0; TTLs.Trial(:,4)];
TrialWiseSniffs(:,10) = odorList(TrialWiseSniffs(:,4),1); % prev. stimulus identity
TrialWiseSniffs(2:end,8) = (TrialWiseSniffs(1:end-1,3)==TrialWiseSniffs(2:end,1)).*(TrialWiseSniffs(1:end-1,3)-TrialWiseSniffs(1:end-1,1)); % prev. sniff duration

%%
% Let's work with an example unit
%whichunit = 48;
unitID = 38;
whichunit = find([SingleUnits.id]==unitID);
window = [0 medianSniff];
sniffwindows = linspace(0,medianSniff,4);
sniffwindows(2,1:3) = sniffwindows(1,2:4);
sniffwindows(:,4) = [];
chunkByPrevOdor = 0;
equalTimeBlocks = 0;
if equalTimeBlocks
    cutOff = [-StimSettings.timing(3)/1000 (StimSettings.timing(3)+StimSettings.timing(4))/1000];
else
    cutOff = [-StimSettings.timing(2)/1000 sum(StimSettings.timing(3:5))/1000];
end

% for every trial, every sniff get a spike count per sniff?
SniffStats = [];
for t = 1:size(TTLs.Trial,1)
    % validsniffs 
    whichsniffs = find( (TrialWiseSniffs(:,4)==t) & (TrialWiseSniffs(:,1)>cutOff(1)) & (TrialWiseSniffs(:,1)<cutOff(2)));
    whichsniffs = TrialWiseSniffs(whichsniffs,:);
    for i = 1:size(whichsniffs,1)
        % spike counts
        thisSniffStart = whichsniffs(i,1) + whichsniffs(i,6); % actual sniff time
        thisSniffEnd   = min((whichsniffs(i,3) + whichsniffs(i,6)), (thisSniffStart + window(2)));
        thisSniffWindow = [thisSniffStart+window(1)  thisSniffEnd];
        % find spikes
        thisSniffSpikes = find((SingleUnits(whichunit).spikes>=thisSniffWindow(1))&(SingleUnits(whichunit).spikes<=thisSniffWindow(2)));
        thisSniffSpikes = SingleUnits(whichunit).spikes(thisSniffSpikes) - thisSniffStart; % actual spiketimes - sniff start
        for s = 1:size(sniffwindows,2)
            spikecounts(s) = numel( find( (thisSniffSpikes>=sniffwindows(1,s)) & (thisSniffSpikes<sniffwindows(2,s)) ) );
        end
        timeFromOdor = thisSniffStart - TTLs.Trial(whichsniffs(i,4),7);
        thisSniffStats = [spikecounts whichsniffs(i,[4 5 7 9 10 8]) whichsniffs(i,3)-whichsniffs(i,1) timeFromOdor];
        SniffStats      = vertcat(SniffStats, thisSniffStats);
    end
end

% SniffStats : Columns
% 1 to 3: SniffCounts in sniffwindows
% 4     : trial index
% 5     : sniff index within a trial, 0 = first sniff after odor onset
% 6     : trial phase : 0 - pre-odor, 1 - odor, 1.5 - second odor pulse, 
%                       2 - post-odor, 2.5 - post second odor pulse
% 7     : Stimulus identity (or mixed if conc. series)
% 8     : Prev Stimulus identity
% 9     : Prev Sniff Duration
% 10    : this Sniff Duration
% 11    : this Sniff time from odor onset


% plot stim wise
stackmode = 1; sniffmax = 6;
figure;
for s = 1:numel(nStim)
    %subplot(numel(nStim),1,s);
    subplot(numel(nStim)/4,numel(nStim)/4,s);
    selectsniffs = find((SniffStats(:,7)==nStim(s))&(abs(SniffStats(:,5))<sniffmax));
    plot(SniffStats(selectsniffs,5),sum(SniffStats(selectsniffs,[1:3]),2),'.','Color',mycolors(s,:));
    H = get(gca,'YLim');
    line([0 0],H,'Color','k','Linewidth',0.5);
    hold on
    plot(SniffStats(selectsniffs,5),sum(SniffStats(selectsniffs,[1:3]),2),'.','Color',mycolors(s,:),'MarkerSize',10);
    if stackmode
        for x = 1:3
            line((sniffmax*2)*x+[0 0],H,'Color','k','Linewidth',0.5);
            selectsniffs = find((SniffStats(:,7)==nStim(s))&(abs(SniffStats(:,5))<6));
            plot((sniffmax*2)*x+ SniffStats(selectsniffs,5),SniffStats(selectsniffs,x),'.','Color',mycolors(s,:),'MarkerSize',10);
        end
    end
end

figure;
for s = 1:numel(nStim)
    %subplot(numel(nStim),1,s);
    subplot(numel(nStim)/4,numel(nStim)/4,s);
%    selectsniffs = find((SniffStats(:,7)==nStim(s))&(abs(SniffStats(:,5))<sniffmax));
    selectsniffs = find((SniffStats(:,7)==nStim(s)));
    %plot(SniffStats(selectsniffs,11),sum(SniffStats(selectsniffs,[1:3]),2),'.','Color',mycolors(s,:));
    plot(SniffStats(selectsniffs,11),sum(SniffStats(selectsniffs,[1:3]),2),'o','MarkerFaceColor','none','MarkerEdgeColor',mycolors(s,:),'MarkerSize',3);
    H = get(gca,'YLim');
    line([0 0],H,'Color','k','Linewidth',0.5);
    hold on    
    % highlight the first sniff
    firstsniffs = find((SniffStats(:,7)==nStim(s))&(SniffStats(:,5)==0));
       plot(SniffStats(firstsniffs,11),sum(SniffStats(firstsniffs,[1:3]),2),'o','MarkerFaceColor',mycolors(s,:),'MarkerEdgeColor','k','MarkerSize',5,'LineWidth',1);
    firstsniffs = find((SniffStats(:,7)==nStim(s))&(SniffStats(:,5)==1));

     plot(SniffStats(firstsniffs,11),sum(SniffStats(firstsniffs,[1:3]),2),'.','Color',mycolors(s,:),'MarkerSize',16);
end

%%

figure;
hold on

sortbyPhase = 1;
phases = [0 1 1.5];
if ~sortbyPhase
    phases = [nan];
end

sortcases = 3;
for subSortCase = 0:1:(sortcases-1)
    %subSortCase = 2; % 0 = just by trial phase, 1 = sort by prev sniff
    %duration, 2 = by sniff index, 3 by previous odor
    sniffsDone = zeros(1,numel(phases));
    for o = 1:numel(nStim)

        for p = 1:numel(phases)
            % select sniffs
            if sortbyPhase
                if chunkByPrevOdor
                    whichsniffs = find( (TrialWiseSniffs(:,10)==nStim(o)) & (TrialWiseSniffs(:,7)==phases(p)));
                else
                    whichsniffs = find( (TrialWiseSniffs(:,9)==nStim(o)) & (TrialWiseSniffs(:,7)==phases(p)));
                end
                subplot(1,numel(phases)*sortcases,p+numel(phases)*subSortCase);
                %            subplot(1,numel(phases),p);
                hold on
            else
                if chunkByPrevOdor
                    whichsniffs = find(TrialWiseSniffs(:,10)==nStim(o));
                else
                    whichsniffs = find(TrialWiseSniffs(:,9)==nStim(o));
                end
            end

            whichsniffs = TrialWiseSniffs(whichsniffs,:);
            % only keep sniffs for which the previous sniff exists
            whichsniffs(whichsniffs(:,8)==0,:) = [];

            if equalTimeBlocks
                if phases(p)==0
                    whichsniffs(whichsniffs(:,1)<cutOff(1),:) = [];
                elseif phases(p)==1.5
                    whichsniffs(whichsniffs(:,1)>cutOff(2),:) = [];
                end
            end

            switch subSortCase
                case 0
                    % sort by phase - irrelevant if sortbyphase == 1
                    whichsniffs = sortrows(whichsniffs,7);
                case 1
                    % sort by phase and duration of previous sniff
                    whichsniffs = sortrows(whichsniffs,[7 8]);
                case 2
                    % sort by trial phase & sniff index 
                    whichsniffs = sortrows(whichsniffs,[7 5]);
                case 3
                    % sort by trial phase & prev. stim
                     whichsniffs = sortrows(whichsniffs,[7 10]);

            end
            % build spikepplot
            SpikesPlot = [];
            for i = 1:size(whichsniffs,1)
                thisSniffStart = whichsniffs(i,1) + whichsniffs(i,6); % actual sniff time
                thisSniffEnd   = min((whichsniffs(i,3) + whichsniffs(i,6)), (thisSniffStart + window(2)));
                thisSniffWindow = [thisSniffStart+window(1)  thisSniffEnd];
                % find spikes
                thisSniffSpikes = find((SingleUnits(whichunit).spikes>=thisSniffWindow(1))&(SingleUnits(whichunit).spikes<=thisSniffWindow(2)));
                thisSniffSpikes = SingleUnits(whichunit).spikes(thisSniffSpikes) - thisSniffStart; % actual spiketimes - sniff start
                SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes i+0*thisSniffSpikes]);
            end
            if ~sortbyPhase
                transitions = [find(diff(whichsniffs(:,7))); i];
                for n = 1:numel(transitions)
                    line(window,sniffsDone(p)+transitions(n)+[0 0],'Color','k');
                end
            else
                line(window,sniffsDone(p)+[0 0],'Color','k');
            end
            plot(SpikesPlot(:,1),(sniffsDone(p))+SpikesPlot(:,2),'.','Color',mycolors(o,:));
            sniffsDone(p) = sniffsDone(p) + i;
        end
    end
    if sortbyPhase
        for p = 1:numel(phases)
            subplot(1,numel(phases)*sortcases,p+numel(phases)*subSortCase);
            %subplot(1,numel(phases),p);
            set(gca,'YLim',[0 sniffsDone(p)]);
            if phases(p)==1
                set(gca,'Box','on','LineWidth',2);
            end
        end
    else
        set(gca,'YLim',[0 sniffsDone]);
    end
end
%%
