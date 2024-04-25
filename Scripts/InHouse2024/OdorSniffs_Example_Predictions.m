%% script to show sniff aligned cells - in ITI
% cells: Q3_20221019_r0_18, Q8_20221204_r0_49, S12_20230731_r0_10,
% S12_20230731_r0_58
warptype = 0;
plotevents = 1;
plotspikes = 1;
psth = 1;
BinOffset = -100; % in ms
shade_exhalation = 1;
plotPSTH = 1;

figure;
for examples = 1:4
    switch examples
        case 1
            SessionTag = 'S12_20230731_r0_processed.mat';   UnitID = 58;    OdorID = 3;
        case 2
            SessionTag = 'S12_20230731_r0_processed.mat';   UnitID = 31;    OdorID = 3;
        case 3
            SessionTag = 'S12_20230731_r0_processed.mat';    UnitID = 44;   OdorID = 3;
        case 4
            SessionTag = 'Q4_20221109_r0_processed.mat';    UnitID = 55;    OdorID = 2;        
    end

MouseName = regexprep(SessionTag,'_(\w+)_processed.mat','');
Paths = WhichComputer();
MySession = fullfile(Paths.ProcessedSessions,MouseName,SessionTag);
PredictedSession = fullfile('/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/sniffPSTHPredictions/', ...
                    MouseName, strrep(SessionTag,'r0_processed','sniffs'));

[TrialAligned, TrialInfo, ...
    ReplayAligned, ReplayInfo, ...
    TuningAligned, TuningInfo, ...
    AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

% get sniff time stamps and info for the sniffs we want to plot
[SelectedSniffs] = SelectSniffs_forKernelFits(TrialAligned, TrialInfo, [1 2 3]);

% some processing
SniffParams = [];
for whichodor = 1:3
    % get the spike raster
    [~, SniffParams_temp] = ...
        SniffAlignedFRs(SelectedSniffs{whichodor}, [], 'onlyparams', 1);
    SniffParams = vertcat(SniffParams, SniffParams_temp);
end
% sniff params
% 1-4: [currsniffstate currsniffloc currinhend currsniffend ...
% 5-6:  currsniffTrialID currsniffIndex ...
% 7:10: prevsniffstate prevsniffloc previnhstart previnhend
% 11:   actual sniff start time in OEPS timebase]
ITISniffs = SniffParams(find(SniffParams(:,1)==-1),:);
tokeep = randperm(size(ITISniffs,1),1000);
ITISniffs = ITISniffs(tokeep,:);
% sort by sniff duration
ITISniffs = sortrows(ITISniffs,4);

OdorSniffs = SniffParams(find(SniffParams(:,1)==OdorID),:);
tokeep = randperm(size(OdorSniffs,1),1000);
OdorSniffs = OdorSniffs(tokeep,:);
% sort by sniff duration
OdorSniffs = sortrows(OdorSniffs,4);

AllSniffs = [ITISniffs; OdorSniffs];

sniffmax(examples)  = max(ITISniffs(:,4)) + 0.2;
thisUnitSpikes = AllUnits.Spikes{UnitID};

% load some stuff from the predictions session
load(PredictedSession,'MyUnits','kernelsout','PSTHbinsize');

% get its predicted PSTH
whichkernel = find(MyUnits==UnitID);
[PSTHOut,kernels{examples},coeffs{examples}] = GetPredictedSniffPSTH(AllSniffs(:,1:10),kernelsout{whichkernel},'binsize', PSTHbinsize);

%% plot the predicted raster
subplot(1,4,examples)
hold on;
dt = 0.010; % in seconds
SpikesPlot = [];
for x = 1:size(AllSniffs,1)   
    ExhalationTimes = [AllSniffs(x,10) 0 AllSniffs(x,3:4)];
    ExhalationTimes = reshape(ExhalationTimes,2,2)';
    if AllSniffs(x,1) < 0
        SniffPlotter(ExhalationTimes', x, Plot_Colors('Odor3'));
    else
        SniffPlotter(ExhalationTimes', x, Plot_Colors('Odor2'));
    end
    
    % get the predicted spike times
    thisSniffSpikes = PechePourPoisson(100*PSTHOut(x,:),dt);
    SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
    %plot(T, x + (0*T), 'k', 'MarkerSize', 2);
end

if ~plotPSTH
    % plot spikes
    plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize', 2);
    %set(gca,'YLim', [0 2200], 'XLim',[-0.05 0.75],'TickDir','out');
    set(gca,'XLim',[-0.05 0.75],'TickDir','out');
else
    imagesc(PSTHOut)
    colormap(brewermap([100],'Blues'))
    set(gca,'YLim',[0 2000],'TickDir','out');
end

end
nSniffs(examples) = size(AllSniffs,1);

set(gcf, 'Position', [281 405 1160 392]);
set(gcf,'Renderer','painters');
if ~plotPSTH
        print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/OdorsniffsExampleRasters_Predictions.eps',...
            '-depsc','-tiff','-r300','-painters');
else
    print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/OdorsniffsExamplePSTHs_Predictions.eps',...
    '-depsc','-tiff','-r300','-painters');
end

if ~plotPSTH
%% plot the kernels
figure;
for examples = 1:4
    subplot(1,4,examples);
    plot(sgolayfilt(kernels{examples}(:,1),1,5),'k');
    hold on
    if examples < 4
        plot(sgolayfilt(kernels{examples}(:,5),1,5),'r');
    else
        plot(sgolayfilt(kernels{examples}(:,4),1,5),'r');
    end
    set(gca,'TickDir','out');   
end

set(gcf, 'Position', [281 405 1160 150]);
set(gcf,'Renderer','painters');

print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/OdorsniffsExampleRasters_kernels.eps',...
    '-depsc','-tiff','-r300','-painters');

end

%%
%handles.(['spikesplot',num2str(i)]) = plot(NaN, NaN, '.k','Markersize', 0.5);
function SniffPlotter(AllTS, rowidx, boxcolor)
%foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.8);
for t = 1:size(AllTS,2)
    TS = AllTS(:,t);
    foo = fill(NaN,NaN,boxcolor,'FaceAlpha',0.7);
    foo.EdgeColor = 'none';
    if ~isempty(TS)
        foo.Vertices = [ ...
            reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
            repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
        foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
    end
end
end
