WhereSession = '/home/priyanka/Desktop/forWDW/S12_20230731_r0_processed.mat';
WhereSession = '/home/priyanka/Desktop/forWDW/Q9_20221116_r0_processed.mat';
WhereSession = '/home/priyanka/Desktop/forWDW/O3_20210927_r0_processed.mat';
%WhereSession = '/home/priyanka/Desktop/forWDW/Q4_20221109_r0_processed.mat';

Paths = WhichComputer();
WhereSession = fullfile(Paths.Wolf.processed,'forWDW','O3_20210927_r0_processed.mat');
WhereSession = fullfile(Paths.Wolf.processed,'forWDW','Q4_20221109_r0_processed.mat');


load(WhereSession);

% loads PassiveOut, TracesOut, SingleUnits
% subfields in *Out: Timestamps, Lever, Motor, Rewards, 
%   Sniffs, SniffsFiltered, SniffsLocationed, SniffsDigitized
%   Trial, Odor, Manifold

%% get sniff timestamps
Sniffidx = find(diff(abs(TracesOut.SniffsDigitized{1})));
Sniffidx = reshape(Sniffidx,2,[])';
Sniffidx(:,1) = Sniffidx(:,1) +  1;
%SniffTS = TracesOut.Timestamps{1}(Sniffidx);

%% separate various sniff types
for n = 1:size(Sniffidx,1)
    % ITI or not
    ITIportion = ...
        numel(find( TracesOut.Manifold{1}(Sniffidx(n,1):Sniffidx(n,2)) )) / ...
        (Sniffidx(n,2) - Sniffidx(n,1) + 1 ); 
    % Odor identity
    StimType = mode( TracesOut.Odor{1}(Sniffidx(n,1):Sniffidx(n,2)) );
    Sniffidx(n,4) = ITIportion;
    Sniffidx(n,5) = StimType;
    Sniffidx(n,6) = mean(TracesOut.SniffsLocationed{1}(Sniffidx(n,1):Sniffidx(n,2))); % odor location
end

%%
thisUnitSpikes = SingleUnits(18).spikes;
window = [-0.1 1];
figure;

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

%% for the passive part
%% get sniff timestamps
Sniffidx = find(diff(abs(PassiveOut.SniffsDigitized{1})));
Sniffidx = reshape(Sniffidx,2,[])';
Sniffidx(:,1) = Sniffidx(:,1) +  1;
%SniffTS = TracesOut.Timestamps{1}(Sniffidx);

%% separate various sniff types
for n = 1:size(Sniffidx,1)
    % ITI or not
    ITIportion = ...
        numel(find( PassiveOut.Manifold{1}(Sniffidx(n,1):Sniffidx(n,2)) )) / ...
        (Sniffidx(n,2) - Sniffidx(n,1) + 1 ); 
    % Odor identity
    StimType = mode( PassiveOut.Odor{1}(Sniffidx(n,1):Sniffidx(n,2)) );
    Sniffidx(n,4) = ITIportion;
    Sniffidx(n,5) = StimType;
    Sniffidx(n,6) = mean(PassiveOut.SniffsLocationed{1}(Sniffidx(n,1):Sniffidx(n,2))); % odor location
    Sniffidx(n,7) = mode( PassiveOut.Trial{1}(Sniffidx(n,1):Sniffidx(n,2)) );
end

%nsniffs = 0;
for snifftype = 1:5
    subplot(2,5,snifftype+5);
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
    whichsniffs = intersect(whichsniffs,find(Sniffidx(:,7)==-4));

    SniffTS = Sniffidx(whichsniffs,:);
    SniffTS(:,1:2) = PassiveOut.Timestamps{1}(SniffTS(:,1:2));
    SniffTS(:,3) = SniffTS(:,2) - SniffTS(:,1);
    SniffTS = sortrows(SniffTS,[6 3],'ascend');
    for x = 1:size(SniffTS,1)
        ts = SniffTS(x,1) + window;
        whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(2)));
        thisSniffSpikes = thisUnitSpikes(whichSpikes) - SniffTS(x,1);
        SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
    end
    if ~isempty(SpikesPlot)
    plot(SpikesPlot(:,1), SpikesPlot(:,2), '.k','Markersize', 0.5);
    %nsniffs = max(nsniffs,x);
    end
end

for snifftype = 1:5
    subplot(2,5,snifftype+5);
    set(gca,'YLim',[0 nsniffs+1]);
end
