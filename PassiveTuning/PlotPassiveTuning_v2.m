function [SpikeCounts,Odors,Locations] = PlotPassiveTuning_v2(AlignedSpikes, TuningTrials, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('savefigures', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('whichunits', [], @(x) isnumeric(x));
params.addParameter('rasters', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('UnitsPerFig', 5, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
savefigs = params.Results.savefigures;
whichUnits = params.Results.whichunits;
plotraster = params.Results.rasters;
plotpsth = params.Results.psth;
nrows = params.Results.UnitsPerFig;

if isempty(whichUnits)
    whichUnits = 1:size(AlignedSpikes,2);
end

Locations = unique(TuningTrials(:,7));
Odors = unique(TuningTrials(:,5));

OdorStart = mean(TuningTrials(find(TuningTrials(:,5)),4) - TuningTrials(find(TuningTrials(:,5)),1));
OdorOff   = mean(TuningTrials(find(TuningTrials(:,5)),6) - TuningTrials(find(TuningTrials(:,5)),1));

ncols = numel(Locations(:,1));
nreps = ceil(size(TuningTrials,1)/numel(Odors)/numel(Locations));
figure;
colortags = {'k', 'r', 'o', 't'};
% Create a matrix of trials

for x = 1:numel(whichUnits) % for every cell
    MyUnit = whichUnits(x);
    
    for MyOdor = 1:numel(Odors)
        thisOdor = Odors(MyOdor);
        for MyLocation = 1:numel(Locations(:,1))
            
            if plotraster
                if mod(x,nrows)
                    whichsubplot = ncols*(mod(x,nrows)-1) + MyLocation;
                    %whichsubplot = mod(MyUnit,20);
                else
                    whichsubplot = ncols*(nrows-1) + MyLocation;
                    %whichsubplot = 20;
                end
                subplot(nrows,ncols,whichsubplot);
                
                if (MyOdor == 1)
                    fill([OdorStart OdorStart OdorOff OdorOff],[0 numel(Odors)*nreps numel(Odors)*nreps 0],[0.8706    0.9216    0.9804],...
                        'EdgeColor','none');
                    hold on;
                end
                
                if MyLocation == 1
                    ylabel(['unit# ',num2str(MyUnit)]);
                end
            end
            
            thisLocation = Locations(MyLocation);
            % get all trial IDs that match this location and this odor
            MyTrials = intersect(find(TuningTrials(:,7)==thisLocation),find(TuningTrials(:,5)==thisOdor));
            
            % get spike counts for each trial
            preodor = []; odor = []; postodor = [];
            for i = 1:numel(MyTrials)
                thisTrialSpikeTimes = AlignedSpikes{MyTrials(i),MyUnit}{1};
%                 % ignore spikes preceding trial start
%                 thisTrialSpikeTimes(thisTrialSpikeTimes<0) = [];
                
                if plotraster
                    % plot the spike raster
                    row_idx = nreps*(MyOdor-1) + i;
                    PlotRaster(thisTrialSpikeTimes,row_idx,Plot_Colors(colortags{MyOdor}));
                end
                % count spikes per stimulus period
                OdorDuration = OdorOff - OdorStart;
                x1 = -OdorDuration; x2 = 0;
                noair(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
                x1 = max(OdorStart-OdorDuration,0); x2 = OdorStart;
                preodor(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
                x1 = OdorStart; x2 = OdorOff;
                odor(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
                x1 = OdorOff; x2 = OdorOff + OdorDuration;
                postodor(i) = (numel(find((thisTrialSpikeTimes>=x1)&(thisTrialSpikeTimes<x2))))/(x2-x1);
            end
            
            SpikeCounts(MyLocation,:,MyUnit,MyOdor) = [mean(noair) mean(preodor) mean(odor) mean(postodor)];           
        end
    end
    
    if plotraster
        % Plot settings
        for MyPlot = whichsubplot-ncols+1:whichsubplot
            subplot(nrows,ncols,MyPlot);
            set(gca,'XTick',[],'YTick',[],'TickDir','out');
        end
        
        
        if mod(x,nrows) == 0
            if savefigs
                saveas(gcf,[MyFileName,'_TuningRaster_',num2str(MyUnit/nrows),'.fig']);
                set(gcf,'renderer','Painters');
                print([MyFileName,'_TuningRaster_',num2str(MyUnit/nrows),'.eps'],'-depsc','-tiff','-r300','-painters');
                close(gcf);
            end
            figure;
        end
    end
end

if savefigs
    saveas(gcf,[MyFileName,'_TuningRaster_',num2str(MyUnit/nrows),'.fig']);
    set(gcf,'renderer','Painters');
    print([MyFileName,'_TuningRaster_',num2str(MyUnit/nrows),'.eps'],'-depsc','-tiff','-r300','-painters');
    close(gcf);
end

% plot tuning curves
figure;
ncols = 1;
for x = 1:numel(whichUnits) % for every cell
    MyUnit = whichUnits(x);
    
    TuningCurve = [];
    for MyOdor = 1:numel(Odors)
       for MyLocation = 1:numel(Locations(:,1))
           TuningCurve(MyOdor,MyLocation) = SpikeCounts(MyLocation,4,MyUnit,MyOdor); 
           %TuningSEM(MyOdor,MyLocation) = SpikeCounts{MyUnit,MyOdor,MyLocation}(2);
       end
    end
    
    if plotpsth
        if mod(x,nrows*ncols)
            subplot(nrows,ncols,mod(x,nrows*ncols));
        else
            subplot(nrows,ncols,nrows*ncols);
        end
        plot(Locations,TuningCurve(1,:),'k');
        hold on
        plot(Locations,TuningCurve(2,:),'color',Plot_Colors('r'));
        plot(Locations,TuningCurve(3,:),'color',Plot_Colors('o'));
        plot(Locations,TuningCurve(4,:),'color',Plot_Colors('t'));
        title(['Unit# ',num2str(MyUnit)]);
        
        
            if mod(x,nrows*ncols) == 0
                if savefigs
                    saveas(gcf,[MyFileName,'_TuningCurve_',num2str(MyUnit/nrows),'.fig']);
                    set(gcf,'renderer','Painters');
                    print([MyFileName,'_TuningCurve_',num2str(MyUnit/nrows),'.eps'],'-depsc','-tiff','-r300','-painters');
                    close(gcf);
                end
                figure;
            end
        
    end
end

if savefigs
    saveas(gcf,[MyFileName,'_TuningCurve_',num2str(MyUnit/nrows),'.fig']);
    set(gcf,'renderer','Painters');
    print([MyFileName,'_TuningCurve_',num2str(MyUnit/nrows),'.eps'],'-depsc','-tiff','-r300','-painters');
    close(gcf);
end

end