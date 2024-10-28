%% script to analyze contamination issues in the O3 dataset

%% load processed sniff and spike data
%SessionName = 'O3_20210927_r0';
SessionName = 'Q4_20221109_r0'; 
MySession = fullfile('/home/priyanka/Desktop/forWDW',[SessionName,'_processed.mat']);
load(MySession); % loads TracesOut PassiveOut SingleUnits

%% Sniff Onsets and Offsets
SniffVector = TracesOut.SniffsDigitized{1};
SniffIdx(:,1) = find(diff([0; SniffVector])>0); % inhalation onset
SniffIdx(:,2) = find(diff([SniffVector; 0])<0); % last inhalation idx

%% lets plot the last sniff in the trial and the next 1-3 sniffs in the ITI
figure;
thisUnitSpikes = SingleUnits(6).spikes;
nSniffs = 10;
%subplot(3,4,n)
for odor = 1:3
    SpikesPlot = []; OdorTransitions = [];
    OdorVector = TracesOut.Odor{1};
    OdorVector(OdorVector~=odor) = 0;
    OdorTransitions(:,1) = find(diff([0; OdorVector])>0); % first index when odor turns ON
    OdorTransitions(:,2) = find(diff([OdorVector; 0])<0); % last index when odor is still ON
    for Transition = 1:size(OdorTransitions,1)-1 % every trial
        % last in trial sniff
        firstsniff = find(SniffIdx(:,1)<OdorTransitions(Transition,2),1,'last');
        nextTrialStart = find(TracesOut.Odor{1}((OdorTransitions(Transition,2)+1):end)>0,1,'first');
        lastsniff = find(SniffIdx(:,1)<OdorTransitions(Transition,2)+nextTrialStart,1,'last');
        for i = firstsniff:lastsniff
            x = i - firstsniff + 1;
            ts = TracesOut.Timestamps{1}(SniffIdx(i,1));
            ts = ts + [-0.25 0 0.75];
            whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(3)));
            thisSniffSpikes = thisUnitSpikes(whichSpikes) - ts(2);
            thisSniffLocation = mean(TracesOut.SniffsLocationed{1}(SniffIdx(i,1):SniffIdx(i,2)));
            thisSniffManifold = mode(TracesOut.Manifold{1}(SniffIdx(i,1):SniffIdx(i,2)));
            SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)...
                Transition*ones(numel(thisSniffSpikes),1) thisSniffLocation*ones(numel(thisSniffSpikes),1) ...
                thisSniffManifold*ones(numel(thisSniffSpikes),1)]);
        end
    end
    for n = 1:nSniffs
        f = find(SpikesPlot(:,2)==n);
        subplot(3,nSniffs,n + (odor-1)*nSniffs); hold on
        mySpikesPlot = SpikesPlot(f,[1 3 4 5]);
        %[~,sortorder] = sortrows(mySpikesPlot,[3 2]);
        mySpikesPlot = sortrows(mySpikesPlot,[4 3 2 1]);
        U = unique(mySpikesPlot(:,2),'stable');
        for x = 1:numel(U)
            mySpikesPlot(find(mySpikesPlot(:,2)==U(x)),5) = x;
        end
        m = find(mySpikesPlot(:,4)==0,1,'last');
        if ~isempty(m)
            mySpikesPlot((m+1):end,5) = mySpikesPlot((m+1):end,5) + 20;
        end
        %plot(SpikesPlot(f,1),SpikesPlot(f,4),'.k');
        f = find(abs(mySpikesPlot(:,3))<10);
        plot(mySpikesPlot(f,1),mySpikesPlot(f,5),'.k');
        f = find(abs(mySpikesPlot(:,3))>=10);
        plot(mySpikesPlot(f,1),mySpikesPlot(f,5),'.r');

        set(gca,'YLim',[0 250]);
    end
end