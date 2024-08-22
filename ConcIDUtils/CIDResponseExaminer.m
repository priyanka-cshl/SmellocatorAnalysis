% load the data
% load('/mnt/grid-hs/mdussauz/CID/Processed/E6/2022-06-10_11-40-39_cid-processed.mat');
% FigPath = '/home/priyanka/Desktop/cid/E6/20220610';

% load('/mnt/grid-hs/mdussauz/CID/Processed/E3/2022-06-14_11-36-14_cid-processed.mat');
% FigPath = '/home/priyanka/Desktop/cid/E3/20220614';

load('/mnt/grid-hs/mdussauz/CID/Processed/E2/2022-06-11_13-57-38_cid-processed.mat');
FigPath = '/home/priyanka/Desktop/cid/E2/20220611';


TTLs.Trial((TTLs.Trial(:,4)==0),4) =  max(TTLs.Trial(:,4)) + 1;
nStim = unique(TTLs.Trial(:,4));
nTypes = unique(TTLs.Trial(:,5));
nreps = max(TTLs.Trial(:,6));
nUnits = size(SingleUnits,2);
savefigs = 1;
mycolors = brewermap(10,'YlOrRd');
mycolors(1:2,:) = [];

%%
for n = 1:nUnits
    if mod(n,10) == 1
        FigureName = ['Units ',num2str(n),'-',num2str(n+9)];
        figure('Name',FigureName);
    end
    row = mod(n,10);
    if ~row
        row = 10;
    end
    thisUnitSpikes = SingleUnits(n).spikes;
    for odor = 1:numel(nStim)
        subplot(10,numel(nStim),((row-1)*numel(nStim) + odor));
        hold on
        for conc = 1:numel(nTypes)
            SpikesPlot = [];
            whichTrials = intersect(find(TTLs.Trial(:,4)==nStim(odor)),find(TTLs.Trial(:,5)==nTypes(conc)));
            for rep = 1:numel(whichTrials)
                ts = TTLs.Trial(whichTrials(rep),[1 2 7 8]); % trial start, stop, odor start, stop
                thistrialspikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<ts(2)));
                thistrialspikes = thisUnitSpikes(thistrialspikes) - ts(3);
                SpikesPlot = vertcat(SpikesPlot, [thistrialspikes (rep + (conc-1)*nreps)*ones(numel(thistrialspikes),1) conc*ones(numel(thistrialspikes),1)]);
            end
            plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(conc*2,:));
        end
        set(gca,'XTick',[],'YTick',[],'XLim',[-6 4],'YLim',[0 0.056]);
        line([0 0],[0 0.056],'color','k');
        line([2 2],[0 0.056],'color','k');
    end

    if ~mod(n,10)
        set(gcf,'Position',[2150 80 1350 900]);
        if savefigs
            saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
            close(gcf);
        end
    end
end

% last figure
if mod(n,10) ~= 0
    set(gcf,'Position',[2150 80 1350 900]);
    if savefigs
        saveas(gcf,fullfile(FigPath,[FigureName,'.png']));
        close(gcf);
    end
end


