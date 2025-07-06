%% script to extract the sniff aligned spiking responses
%  of a given neuron 

%% Step 1: Data Paths
myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31';

%% Get the sniff timestamps
%  Load Sniffs from the T3 session and make a vector with binarized sniffs
load(fullfile(myKsDir,'quickprocesssniffs.mat'));
load(fullfile(myKsDir,'SessionDetails.mat'));
startTS = 0;
for i = 1:numel(Files.Samples)
    endTS = Files.Samples(i)/30000 + startTS; % in seconds
    whichsniffs = find( (AllSniffs(:,1)>=startTS) & (AllSniffs(:,1)<endTS) );
    AllSniffs(whichsniffs,8) = i; % session phase
    startTS = endTS;
end

% group sniffs by odors
% there are no air OFF sniffs here
AllSniffs(:,4) = 1; % air is always on
[ParsedSniffs, StimulusList] = ParseSniffsByStimuli(AllSniffs, 'SortBy', 3);
Sniffs2Use{1} = ParsedSniffs{2}; % Air sniffs
% Sniffs2Use{1}(:,8) = 1;
Sniffs2Use{1}(11000:end,:) = [];
Sniffs2Use{1}(5000:8500,:) = [];

%% Get the PSTHs at the same resolution
%  Load spikes
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

%selectedUnits = [19 24 29 30 54 56 57 101 108 133]; % 134 167
selectedUnits = [19 54 24 56 57 29 30 108 101 133];

%%
savefigs = 1;
nRows = 6*2; nCols = 10; 

for thisunit = 1:numel(selectedUnits);
    
    if ~rem(thisunit-1,nCols);
        figure;
    end

    % actual spikes
    k = selectedUnits(thisunit);
    thisUnitSpikes = SingleUnits(k).spikes;
    % sort by duration only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},3,'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    whichplots = (1:10:(nRows*nCols*0.5)) + rem(thisunit-1,nCols);
    subplot(nRows+1,nCols,whichplots);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    
    % sort by session phase and duration only
    Sniffs2Use{1} = sortrows(Sniffs2Use{1},[8 3],'ascend');
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    whichplots = whichplots + (nRows*nCols*0.5);
    subplot(nRows+1,nCols,whichplots);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);

    if ~rem(thisunit,nCols)
        set(gcf,'Position',[2015 535 1748 323]);
        if savefigs
            set(gcf,'Color','w');
            set(gcf,'renderer','Painters');
            print(['/mnt/data/simulation/ObservedAirResponse',num2str(thisunit),'.eps'],'-depsc','-tiff','-r300','-painters');
            close(gcf);
        end
    end
    
    
end
