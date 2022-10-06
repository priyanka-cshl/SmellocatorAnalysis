function [SpikeRasterOut,SpikeCounts,Odors,Locations] = ForPassiveTuningDecoder(AlignedSpikes, TuningTrials, window)

whichUnits = 1:size(AlignedSpikes,2);

Locations = unique(TuningTrials(:,7));
Odors = unique(TuningTrials(:,5));

OdorStart = mean(TuningTrials(find(TuningTrials(:,5)),4) - TuningTrials(find(TuningTrials(:,5)),1));
OdorOff   = mean(TuningTrials(find(TuningTrials(:,5)),6) - TuningTrials(find(TuningTrials(:,5)),1));

% Dimensions of output matrix
% Neurons x odor x location x timebins x repeats

for x = 1:numel(whichUnits) % for every cell
    MyUnit = whichUnits(x);
    
    for MyOdor = 1:numel(Odors)
        thisOdor = Odors(MyOdor);
        
        for MyLocation = 1:numel(Locations(:,1))
            thisLocation = Locations(MyLocation);
            % get all trial IDs that match this location and this odor
            MyTrials = intersect(find(TuningTrials(:,7)==thisLocation),find(TuningTrials(:,5)==thisOdor));
            
            % get spike counts for each trial
            preodor = []; odor = []; postodor = [];
            for i = 1:numel(MyTrials)
                thisTrialSpikeTimes = AlignedSpikes{MyTrials(i),MyUnit}{1};
                
                % count spikes per stimulus epoch
                OdorDuration = OdorOff - OdorStart;
                x1 = -OdorDuration; x2 = 0;
                noair(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
                x1 = max(OdorStart-OdorDuration,0); x2 = OdorStart;
                preodor(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
                x1 = OdorStart; x2 = OdorOff;
                odor(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
                x1 = OdorOff; x2 = OdorOff + OdorDuration;
                postodor(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
                
                % make a raster
                myRaster = zeros(1,diff(window)*1000);
                % align to window start
                thisTrialSpikeTimes = thisTrialSpikeTimes - window(1);
                % convert spike times to milliseconds and floor values
                thisTrialSpikeTimes = ceil(1000*thisTrialSpikeTimes);
                % remove NaNs
                thisTrialSpikeTimes(isnan(thisTrialSpikeTimes)) = [];
                % Make raster
                [C,~,ic] = unique(thisTrialSpikeTimes);
                bin_counts = accumarray(ic,1);
                if ~isempty(C)
                    % ignore any time bins less than 1 sec before trial start
                    bin_counts((C<=0),:) = [];
                    C(C<=0) = [];
                    myRaster(1,C) = bin_counts; %#ok<AGROW>
                end
                
                SpikeRasterOut(x,MyOdor,MyLocation,1:numel(myRaster),i) = myRaster;

            end
            
            SpikeCounts(:,MyLocation,x,MyOdor) = [mean(noair) mean(preodor) mean(odor) mean(postodor)];           
        end
    end
end

end