function [] = PlotLocationTuningRaster(RasterIn, WhichUnit, WhichLocation)
colortags = {'k', 'r', 'o', 't'};
rownum = 0;
for odor = 1:size(RasterIn,2)
    for rep = 1:size(RasterIn,5)
        rownum = rownum + 1;
        spiketimes = find(squeeze(RasterIn(WhichUnit,odor,WhichLocation,:,rep)))/1000; 
        PlotRaster(spiketimes,rownum,Plot_Colors(colortags{odor})); 
    end
end
end