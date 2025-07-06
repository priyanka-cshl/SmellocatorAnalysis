function [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(GroupedSniffs, thisUnitSpikes, yoffset)
if nargin < 3
    yoffset = 0;
end

window = [-0.1 1];

maxsniffs = 0;
for snifftype = 1:size(GroupedSniffs,2)
    SpikesPlot = [];
    SniffTS = GroupedSniffs{snifftype};

    % plot spikes
    for x = 1:size(SniffTS,1)
        ts = SniffTS(x,1) + window;
        whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(2)));
        thisSniffSpikes = thisUnitSpikes(whichSpikes) - SniffTS(x,1);
        y = x + yoffset;
%         if SniffTS(x,8) == 1
            SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes y*ones(numel(thisSniffSpikes),1)]);
%         else
%             SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes -(y+100*(SniffTS(x,8)-1))*ones(numel(thisSniffSpikes),1)]);
%         end
    end

    %plot(SpikesPlot(:,1), SpikesPlot(:,2), '.k','Markersize', 0.5);
    SpikeRaster{snifftype} = SpikesPlot;
    maxsniffs = max(maxsniffs,x);
end

end