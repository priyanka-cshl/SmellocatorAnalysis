%% Settings
PSTHbinsize = 20; 

%% Load Sniffs from the T3 session and make a vector with binarized sniffs
load('/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31/quickprocesssniffs.mat');
SniffTrace(:,1) = RespirationData(:,1); % Timestamps in 500hz base
SniffTrace(:,2) = 0;

for i = 1:size(AllSniffs,1) % every sniff timestamp
    idx = AllSniffs(i,11:12);
    SniffTrace(idx(1):idx(2),2) = 1;
end

% figure; 
% plot(RespirationData(:,1),RespirationData(:,3));
% hold on
% plot(SniffTrace(:,1),SniffTrace(:,2));

% downsample to 20ms binsize
downsample  = PSTHbinsize/(1000/500);
TS_temp     = SniffTrace(:,1)';
TS_down     = TS_temp(1):PSTHbinsize/1000:TS_temp(end);

SniffTraceDS(:,1) = TS_down;
SniffTraceDS(:,2) = interp1(TS_temp,SniffTrace(:,2)',TS_down,'nearest');

%% some extra sniff processing for later
load('/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31/SessionDetails.mat');
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
Sniffs2Use{1}(:,8) = 1;

%% Load Kernels
load('/mnt/data/MTKernels.mat');

%% figures
savefigs = 1;
nRows = 6; nCols = 6;

for k = 1:size(KernelsSmooth,1) % every kernel
    ts = (1:41)*0.05; % 20hz resolution
    ts_new = 0.02:0.02:ts(end);
    KernelDS = interp1(ts,KernelsSmooth(k,:),ts_new);
    KernelDS(:,find(isnan(KernelDS))) = [];

    % convolve with the sniff trace
    PSTHout = conv(SniffTraceDS(:,2)',KernelDS,'full');
    PSTHout(PSTHout<0) = 0;
    %PSTHout(:,1:length(KernelDS)) = [];
    thisUnitSpikes = PechePourPoisson(100*PSTHout,PSTHbinsize/1000);
    [SpikeRaster, maxsniffs] = GetSniffLockedSpikes(Sniffs2Use, thisUnitSpikes);
    
    if ~rem(k-1,nCols)
        FigureName = ['Units ',num2str(k),' - ',num2str(k+nCols-1)]; % one figure per cell
        figure('Name',FigureName);
    end
    
    whichcol = rem(k-1,nCols) + 1;
    subplot(nRows,nCols,whichcol);
    plot(KernelDS);
    
    whichplots = (7:nCols:(nRows*nCols)) + (whichcol-1);
    subplot(nRows,nCols,whichplots);
    plot(SpikeRaster{1}(:,1), SpikeRaster{1}(:,2), '.k','Markersize', 0.5);
    
    if ~rem(k,nCols)
        set(gcf,'Position',[597   234   966   704]);
        if savefigs
            set(gcf,'Color','w');
            set(gcf,'renderer','Painters');
            print(['/mnt/data/simulation/snifftuning',num2str(k),'.eps'],'-depsc','-tiff','-r300','-painters');
            close(gcf);
        end
    end
end


