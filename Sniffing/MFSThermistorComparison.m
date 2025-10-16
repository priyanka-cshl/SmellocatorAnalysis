
WhereSession = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31/quickprocesssniffs.mat';
noMFSafter = 6304.7;

WhereSession = '/mnt/data/Sorted/T2/_2025-05-21_09-18-56_2025-05-21_11-11-12_2025-05-21_11-23-51/quickprocesssniffs.mat';
noMFSafter = 6403;

load(WhereSession,"SniffCoords","SniffProps");
noMFSafter = find(SniffCoords(:,1)>noMFSafter,1,"first");

SniffCoords(noMFSafter:end,:) = [];

SniffCoords(:,16) = SniffCoords(:,3)-SniffCoords(:,1);
SniffCoords = sortrows(SniffCoords,16);
ThermPeaks = SniffCoords(:,1);
MFS2ThermPeaks = SniffCoords(:,6);
MFSCrossings = SniffCoords(:,11);

SniffDurations = SniffCoords(:,16);

%%
maxY = 0.15;

%% thermistor peaks vs mfs peaks
subplot(2,3,1)
scatter(SniffDurations,ThermPeaks-MFS2ThermPeaks,'.');
set(gca,'YLim', [-0.05 maxY]);

% bin sniffs
binnedVals = [];
f1 = 1;
for x = 0.1:0.1:0.99
    f2 = find(SniffCoords(:,16)<=x,1,'last');
    binnedVals = vertcat(binnedVals, ...
        [x median(ThermPeaks(f1:f2)-MFS2ThermPeaks(f1:f2)) std(ThermPeaks(f1:f2)-MFS2ThermPeaks(f1:f2))] );
    f1 = f2+1;
end
subplot(2,3,4)
errorbar(binnedVals(:,1),binnedVals(:,2),binnedVals(:,3))

%% thermistor peaks vs mfs zero crossings
subplot(2,3,2)
scatter(SniffDurations,ThermPeaks-MFSCrossings,'.');
set(gca,'YLim', [-0.05 maxY]);

% bin sniffs
binnedVals = [];
f1 = 1;
for x = 0.1:0.1:0.99
    f2 = find(SniffCoords(:,16)<=x,1,'last');
    binnedVals = vertcat(binnedVals, ...
        [x median(ThermPeaks(f1:f2)-MFSCrossings(f1:f2)) std(ThermPeaks(f1:f2)-MFSCrossings(f1:f2))] );
    f1 = f2+1;
end
subplot(2,3,5)
errorbar(binnedVals(:,1),binnedVals(:,2),binnedVals(:,3))

%% mfs peaks vs mfs zero crossings
subplot(2,3,3)
scatter(SniffDurations,MFS2ThermPeaks-MFSCrossings,'.');
set(gca,'YLim', [-0.05 maxY]);

% bin sniffs
binnedVals = [];
f1 = 1;
for x = 0.1:0.1:0.99
    f2 = find(SniffCoords(:,16)<=x,1,'last');
    binnedVals = vertcat(binnedVals, ...
        [x median(MFS2ThermPeaks(f1:f2)-MFSCrossings(f1:f2)) std(MFS2ThermPeaks(f1:f2)-MFSCrossings(f1:f2))] );
    f1 = f2+1;
end
subplot(2,3,6)
errorbar(binnedVals(:,1),binnedVals(:,2),binnedVals(:,3))
