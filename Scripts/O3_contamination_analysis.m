%% script to analyze contamination issues in the O3 dataset

%% load processed sniff and spike data
SessionName = 'O3_20210927_r0';
MySession = fullfile('/home/priyanka/Desktop/forWDW',[SessionName,'_processed.mat']);
load(MySession); % loads TracesOut PassiveOut SingleUnits

%% Sniff Onsets and Offsets
SniffVector = TracesOut.SniffsDigitized{1};
SniffIdx(:,1) = find(diff([0; SniffVector])>0); % inhalation onset
SniffIdx(:,2) = find(diff([SniffVector; 0])<0); % last inhalation idx

%% lets plot the last sniff in the trial and the next 1-3 sniffs in the ITI
figure;
thisUnitSpikes = SingleUnits(72).spikes;
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
            SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1) Transition*ones(numel(thisSniffSpikes),1)]);
        end
    end
    for n = 1:4
        f = find(SpikesPlot(:,2)==n);
        subplot(3,4,n + (odor-1)*4); 
        plot(SpikesPlot(f,1),SpikesPlot(f,3),'.k');
    end
end