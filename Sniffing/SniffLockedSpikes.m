function [] = SniffLockedSpikes(MySession, whichUnit)

%% load processed sniff and spike data
load(MySession); % loads TracesOut PassiveOut SingleUnits

%% Sniff Onsets and Offsets
Sniffidx = find(diff(abs(TracesOut.SniffsDigitized{1})));
Sniffidx = reshape(Sniffidx,2,[])';
Sniffidx(:,1) = Sniffidx(:,1) +  1;
Sniffidx(1:end-1,3) = Sniffidx(2:end,1); % next inhalation index
Sniffidx(end,:) = [];

%% separate various sniff types
for n = 1:size(Sniffidx,1)

    % ITI or not
    ITIportion = ...
        numel(find( TracesOut.Manifold{1}(Sniffidx(n,1):Sniffidx(n,2)) )) / ...
        (Sniffidx(n,2) - Sniffidx(n,1) + 1 ); % fraction of time when manifold was on during inhalation
    % Odor identity
    StimType = mode( TracesOut.Odor{1}(Sniffidx(n,1):Sniffidx(n,2)) );
    Sniffidx(n,4) = ITIportion;
    Sniffidx(n,5) = StimType;
    Sniffidx(n,6) = mean(TracesOut.SniffsLocationed{1}(Sniffidx(n,1):Sniffidx(n,2))); % odor location
    Sniffidx(n,7) = mode( TracesOut.Trial{1}(Sniffidx(n,1):Sniffidx(n,2)) ); % perturbed or not 
    
end

%% plotting sniff-locked spikes
%whichUnit = 2;
thisUnitSpikes = SingleUnits(whichUnit).spikes;
window = [-0.1 1];
figure;

%% split by stimulus conditions 
nsniffs = 0;
for snifftype = 1:5
    subplot(2,5,snifftype);
    SpikesPlot = [];

    switch snifftype
        case 1 % ITI
            whichsniffs = find(Sniffidx(:,4)==0);
        case 2
            whichsniffs = intersect(find(Sniffidx(:,4)==1), find(Sniffidx(:,5)==0));
        case 3
            whichsniffs = intersect(find(Sniffidx(:,4)==1), find(Sniffidx(:,5)==1));
        case 4
            whichsniffs = intersect(find(Sniffidx(:,4)==1), find(Sniffidx(:,5)==2));
        case 5
            whichsniffs = intersect(find(Sniffidx(:,4)==1), find(Sniffidx(:,5)==3));
    end

    SniffTS = Sniffidx(whichsniffs,:);
    SniffTS(:,1:2) = TracesOut.Timestamps{1}(SniffTS(:,1:2));
    SniffTS(:,3) = SniffTS(:,2) - SniffTS(:,1);
    SniffTS = sortrows(SniffTS,3,'ascend');
    for x = 1:size(SniffTS,1)
        ts = SniffTS(x,1) + window;
        whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(2)));
        thisSniffSpikes = thisUnitSpikes(whichSpikes) - SniffTS(x,1);
        SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
    end
    plot(SpikesPlot(:,1), SpikesPlot(:,2), '.k','Markersize', 0.5);
    nsniffs = max(nsniffs,x);
end

for snifftype = 1:5
    subplot(2,5,snifftype);
    set(gca,'YLim',[0 nsniffs+1]);
end