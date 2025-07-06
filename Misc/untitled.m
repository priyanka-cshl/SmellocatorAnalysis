myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
figure, plot(RespirationData(:,1), RespirationData(:,3));
SingleUnits = GetSingleUnits(myKsDir, 3);
selectedUnits = [19 54 24 56 57 29 30 108 101 133];
hold on; 
for i = 1:10
    PlotRaster(SingleUnits(selectedUnits(i)).spikes,-i/2,'k',0.25);

    %plot(SingleUnits(selectedUnits(i)).spikes, -i/4, '.k','Markersize', 0.5); 
end