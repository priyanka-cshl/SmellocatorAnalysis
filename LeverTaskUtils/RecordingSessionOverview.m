function [] = RecordingSessionOverview(SingleUnits)
    figure; hold on
    
    for i = 1:size(SingleUnits,2)
        foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id];
    end
    [~,SortedByTetrodes] = sort(foo(:,1));
    ClusterIds = foo(SortedByTetrodes,2);
    whichtetrode = SingleUnits(SortedByTetrodes(1)).tetrode;
    MyColors(1,:) = Plot_Colors('r');
    MyColors(2,:) = Plot_Colors('k');
    for i = 1:1:size(SingleUnits,2)
        i
        whichunit = SortedByTetrodes(i);
        thistetrode = SingleUnits(whichunit).tetrode;
        if thistetrode~=whichtetrode
            MyColors = circshift(MyColors,1);
            whichtetrode = thistetrode;
        end
        PlotRaster(SingleUnits(whichunit).spikes,i,MyColors(1,:),0.8);
    end
end