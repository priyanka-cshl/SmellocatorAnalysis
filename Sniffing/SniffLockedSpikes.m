function [] = SniffLockedSpikes(MySession, whichUnit)

%% load processed sniff and spike data
[AllSniffs, ColumnInfo, SingleUnits] = GetAllSniffs(MySession);

%% plotting sniff-locked spikes
%whichUnit = 2;
thisUnitSpikes = SingleUnits(whichUnit).spikes;
figure;

%% split by stimulus conditions 
GroupedSniffs = ParseSniffsByType(AllSniffs);

[SpikeRaster, maxsniffs] = GetSniffLockedSpikes(GroupedSniffs, thisUnitSpikes);

for snifftype = 1:5
    subplot(1,5,snifftype);
    SpikesPlot = SpikeRaster{snifftype};
    SpikesPlot(:,2) = abs(SpikesPlot(:,2));
    plot(SpikesPlot(:,1), SpikesPlot(:,2), '.k','Markersize', 0.5);
    set(gca,'YLim',[0 maxsniffs+1]);
end