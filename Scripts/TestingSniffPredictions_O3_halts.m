%% load the sniff prediction kernels, fits etc.
load('/home/priyanka/Desktop/sniffPSTHPredictions/O3/O3_20210927_sniffs.mat');
% loads MyUnits, PSTHbinsize, SessionName, Sniffs, kernelsout

nUnits = size(kernelsout,1);

plotPSTH = 0;
plotRasters = 0;
plothaltRasters = 1;
plotResiduals = 0;
BinOffset = -100; % in ms
print2pdf = 1;

for whichUnit = 1:nUnits %5; % original unit name
    
    % for any unit - there are two sets of sniffs
    unitID = find(MyUnits==whichUnit); % index in the fitted dataset
    
    Ob_PSTHs = Sniffs.PSTHs{unitID};
    Tr_PSTHs = zeros(size(Ob_PSTHs,1),size(Ob_PSTHs,2)-1); % for predictions using the training kernel
    Cv_PSTHs = Tr_PSTHs; % for predictions using the crossvalidated kernel
    for sniffset = 1:2
        whichsniffs = Sniffs.IDs{unitID}{sniffset}; % which sniffs were used to get the kernel
        sniffparams = Sniffs.Params{unitID}(whichsniffs,:);
        
        kernelID = sniffset; % training kernel
        [PSTHOut] = GetPredictedSniffPSTH(sniffparams,kernelsout{unitID,kernelID},'binsize', PSTHbinsize);
        Tr_PSTHs(whichsniffs,1:size(PSTHOut,2)) = PSTHOut;
        
        kernelID = 1 + mod(sniffset,2); % crossvalidation kernel
        [PSTHOut] = GetPredictedSniffPSTH(sniffparams,kernelsout{unitID,kernelID},'binsize', PSTHbinsize);
        Cv_PSTHs(whichsniffs,1:size(PSTHOut,2)) = PSTHOut;
    end
    
    % sort sniffs by odor ID and then by duration?
    [AllSniffs,sortorder] = sortrows(Sniffs.Params{unitID},[1 4]); % by odor state and by sniff duration
    Ob_PSTHs = Ob_PSTHs(sortorder,:);
    Tr_PSTHs = Tr_PSTHs(sortorder,:);
    Cv_PSTHs = Cv_PSTHs(sortorder,:);
    
    % also rectify FRs
    Tr_PSTHs(Tr_PSTHs<0) = 0;
    Cv_PSTHs(Cv_PSTHs<0) = 0;
    
    %% Calculate fit quality
    Res_training = (Ob_PSTHs(:,2:end) - Tr_PSTHs).^2;
    Res_crossval = (Ob_PSTHs(:,2:end) - Cv_PSTHs).^2;
    for n = 1:size(Ob_PSTHs,1) % every sniff
        % take the invalid points to NaN
        bin_idx = 1 + Ob_PSTHs(n,1);
        if numel(Res_training(n,:)) >= bin_idx
            Res_training(n,bin_idx:end) = NaN;
        end
        if numel(Res_crossval(n,:)) >= bin_idx
            Res_crossval(n,bin_idx:end) = NaN;
        end
    end
    
    if plotResiduals
        % plot residuals
        figure;
        histogram(mean(Res_training,2,'omitnan'),'DisplayStyle','stairs','EdgeColor','k')
        hold on; histogram(mean(Res_crossval,2,'omitnan'),'DisplayStyle','stairs','EdgeColor','r')
        set(gca,'XLim',[-0.05 1], 'YLim', [0 5000]);
    end
    
    % keep track of residuals
    Residuals(unitID,1) = mean(mean(Res_training,2,'omitnan'));
    Residuals(unitID,2) = mean(mean(Res_crossval,2,'omitnan'));
    
    %%
    % load the actual spikes
    MouseName = SessionName(1:regexp(SessionName,'_','once')-1);
    SessionTag = [SessionName,'_processed.mat'];
    Paths = WhichComputer();
    MySession = fullfile(Paths.ProcessedSessions,MouseName,SessionTag);
    
    [TrialAligned, TrialInfo, ...
        ReplayAligned, ReplayInfo, ...
        TuningAligned, TuningInfo, ...
        AllUnits] = ...
        PreprocessSpikesAndSniffs(MySession);
    thisUnitSpikes = AllUnits.Spikes{unitID};

    %%
    
    % plot the observed raster
    % plot the spike rasters
    
    % need the actual OEPS time of the sniff for raw spike raster plotting
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
    SniffParams = sortrows(SniffParams,[1 4]);
    if ~any(any(SniffParams(:,1:10) - AllSniffs))
        AllSniffs(:,11) = SniffParams(:,11);
    end
    
    if plotRasters
        
        figure;
        % plot the observed and predicted rasters
        for type = 1:3
            switch type
                case 1
                    subplot(1,3,1);
                    hold on
                    SniffAlignedSpikesPlotter(AllSniffs, thisUnitSpikes, [], PSTHbinsize/1000, BinOffset);
                case 2
                    % training fits
                    subplot(1,3,2);
                    hold on
                    SniffAlignedSpikesPlotter(AllSniffs, [], Tr_PSTHs, PSTHbinsize/1000, BinOffset);
                case 3
                    % crossvalidated fits
                    subplot(1,3,3);
                    hold on
                    SniffAlignedSpikesPlotter(AllSniffs, [], Cv_PSTHs, PSTHbinsize/1000, BinOffset);
            end
        end
        
    end
    
    
    %% get perturbation sniffs
    [AllPerturbationSniffs] = SelectSniffs_forPerturbationFits(TrialAligned, TrialInfo);
    % get the oberved PSTHs
    [~, PerturbationSniffParams] = SniffAlignedFRs(AllPerturbationSniffs, [], 'onlyparams', 1);
    [PerturbationPSTHs, ~] = SniffAlignedFRs(AllPerturbationSniffs, thisUnitSpikes, ...
                    'psthtype', 'binned', 'binsize', PSTHbinsize);
    % sort by sniff duration and location
    [PerturbationSniffParams,p_sort] = sortrows(PerturbationSniffParams,[1 4]);
    PerturbationPSTHs = PerturbationPSTHs(p_sort,:);
    
    % get the predicted PSTHs
    [PerturbationPSTH_1] = GetPredictedSniffPSTH(PerturbationSniffParams(:,1:10),kernelsout{unitID,1},'binsize', PSTHbinsize);
    [PerturbationPSTH_2] = GetPredictedSniffPSTH(PerturbationSniffParams(:,1:10),kernelsout{unitID,2},'binsize', PSTHbinsize);
    
    % rectify
    PerturbationPSTH_1(PerturbationPSTH_1<0) = 0;
    PerturbationPSTH_2(PerturbationPSTH_2<0) = 0;
    
    % do the same above steps for control sniffs
    % get a set of matched sniffs (controls)
    % similar locations
    matchedsniffs = find((SniffParams(:,1)==PerturbationSniffParams(1,1))&(abs(round(SniffParams(:,2)-30))<5));
    LocationMatchedSniffs = SniffParams(matchedsniffs,:);
    % sort by location and duration
    LocationMatchedSniffs = sortrows(LocationMatchedSniffs,[1 4]);
    
    % get the predicted PSTHs
    [LocationMatchedPSTH_1] = GetPredictedSniffPSTH(LocationMatchedSniffs(:,1:10),kernelsout{unitID,1},'binsize', PSTHbinsize);
    [LocationMatchedPSTH_2] = GetPredictedSniffPSTH(LocationMatchedSniffs(:,1:10),kernelsout{unitID,2},'binsize', PSTHbinsize);
    
    % plot the rasters
    if plothaltRasters
        %figure;
        foo = ['Unit ',num2str(unitID)];
        figure('Name',foo);
        
        % observed
        for type = 1:6
            switch type
                case 1 % halt sniffs
                    subplot(4,2,1);
                    hold on
                    haltPSTH = SniffAlignedSpikesPlotter(PerturbationSniffParams, thisUnitSpikes, [], PSTHbinsize/1000, BinOffset);
                    subplot(4,2,7); 
                    hold on
                    plot(haltPSTH(:,1),haltPSTH(:,2),'k');
                    set(gca,'XLim',[-0.05 0.75],'TickDir','out'); 
                    myYLim = get(gca,'YLim');
                    myYLim(2) = myYLim(2) + 10;
                    set(gca,'YLim', myYLim);
                case 2
                    subplot(4,2,3);
                    hold on
                    SniffAlignedSpikesPlotter(PerturbationSniffParams, [], PerturbationPSTH_1, PSTHbinsize/1000, BinOffset);
                    [meanPSTH] = AverageSniffPSTH(PerturbationPSTH_1,PerturbationSniffParams(:,4), PSTHbinsize);
                    subplot(4,2,7); 
                    plot(meanPSTH(:,1),meanPSTH(:,2),'r');
                    set(gca,'XLim',[-0.05 0.75],'YLim',myYLim,'TickDir','out');   
                case 3
                    subplot(4,2,5);
                    hold on
                    SniffAlignedSpikesPlotter(PerturbationSniffParams, [], PerturbationPSTH_2, PSTHbinsize/1000, BinOffset);
                case 4
                    subplot(4,2,2);
                    hold on
                    ctrlPSTH = SniffAlignedSpikesPlotter(LocationMatchedSniffs, thisUnitSpikes, [], PSTHbinsize/1000, BinOffset);
                    subplot(4,2,8); 
                    hold on
                    plot(ctrlPSTH(:,1),ctrlPSTH(:,2),'k');
                    set(gca,'XLim',[-0.05 0.75],'YLim',myYLim,'TickDir','out');
                case 5
                    subplot(4,2,4);
                    hold on
                    SniffAlignedSpikesPlotter(LocationMatchedSniffs, [], LocationMatchedPSTH_1, PSTHbinsize/1000, BinOffset);
                    [meanPSTH] = AverageSniffPSTH(LocationMatchedPSTH_1,LocationMatchedSniffs(:,4), PSTHbinsize);
                    subplot(4,2,8); 
                    plot(meanPSTH(:,1),meanPSTH(:,2),'r');
                    set(gca,'XLim',[-0.05 0.75],'YLim',myYLim,'TickDir','out');
                case 6
                    subplot(4,2,6);
                    hold on
                    SniffAlignedSpikesPlotter(LocationMatchedSniffs, [], LocationMatchedPSTH_2, PSTHbinsize/1000, BinOffset);
            end
        end
        
        if print2pdf
            FigPath = ['/home/priyanka/Desktop/haltPredictions/', SessionName(1:regexp(SessionName,'_','once')-1)];
            if ~exist(FigPath,'dir')
                mkdir(FigPath);
                fileattrib(FigPath, '+w','a');
            end
            saveas(gcf,fullfile(FigPath,[foo,'.png']),'png');
            close(gcf);
        end
        
    end
    
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

function [myPSTH] = SniffAlignedSpikesPlotter(AllSniffParams, thisUnitSpikes, PSTHIn, dt, BinOffset)

SpikesPlot = [];
for x = 1:size(AllSniffParams,1)
% Plot inhalation periods - adjust times if needed
%     InhalationTimes = [AllSniffParams(x,[9 10]) 0 AllSniffParams(x,3)];
%     InhalationTimes = reshape(InhalationTimes,2,2)';

% Plot exhalation periods - adjust times if needed
    ExhalationTimes = [AllSniffParams(x,10) 0 AllSniffParams(x,3:4)];
    ExhalationTimes = reshape(ExhalationTimes,2,2)';
    
    if AllSniffParams(x,1) < 0
        SniffPlotter(ExhalationTimes', x, Plot_Colors('Odor3'));
    else
        SniffPlotter(ExhalationTimes', x, Plot_Colors('Odor2'));
    end
    
    if ~isempty(thisUnitSpikes)
        ts = [AllSniffParams(x,9) 0 AllSniffParams(x,4) AllSniffParams(x,3)] + AllSniffParams(x,11);
        % [prevsniffstart thisniffstart thissniffend thisinhend]
        whichSpikes = intersect(find(thisUnitSpikes>=ts(2)), find(thisUnitSpikes<=ts(3)));
        
        % align to inhalation start
        thisSniffSpikes = thisUnitSpikes(whichSpikes) - ts(2);
    else
        % get the predicted spike times from the PSTH
        thisSniffSpikes = PechePourPoisson(100*PSTHIn(x,:),dt);
    end
    
    SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
end

% plot spikes
plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize', 2);
%set(gca,'YLim', [0 2200], 'XLim',[-0.05 0.75],'TickDir','out');
set(gca,'XLim',[-0.05 0.75],'TickDir','out');

% make the PSTH
myPSTH = MakeSniffTriggeredPSTH(SpikesPlot(:,1),AllSniffParams(:,4),...
                 BinOffset,'downsample',500,'kernelsize',20);
ts = (1:length(myPSTH))*0.002;
ts = ts + (BinOffset/1000);
myPSTH = [ts' myPSTH];
end

function [meanPSTH] = AverageSniffPSTH(allPSTH,PSTHdurations, PSTHbinsize)
PSTHlengths = floor(PSTHdurations*(1000/PSTHbinsize));
for i = 1:numel(PSTHlengths)
    if PSTHlengths(i)<=size(allPSTH,2)
        allPSTH(i,(PSTHlengths(i)+1):end) = NaN;
    end
end
meanPSTH = mean(allPSTH,1,'omitnan');
ts = (1:length(meanPSTH))*PSTHbinsize/1000;
meanPSTH = [ts' 100*meanPSTH'];
end

