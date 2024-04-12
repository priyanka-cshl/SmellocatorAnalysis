%% script to extract the sniff aligned closed-loop spiking responses
%  of a given neuron and fit the ITI, Air and Odor kernels

%% Step 1: Get the spiking data and sniff parameters
%SessionName = 'Q4_20221109_r0';
%SessionName = 'S12_20230731_r0';
SessionName = 'Q9_20221116_r0';
MySession = fullfile('/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/', ...
                        SessionName(1:regexp(SessionName,'_','once')-1), ...
                        [SessionName,'_processed.mat']);
FigPath = ['/home/priyanka/Desktop/sniffPSTHPredictions/', SessionName(1:regexp(SessionName,'_','once')-1)];
if ~exist(FigPath,'dir')
    mkdir(FigPath);
    fileattrib(FigPath, '+w','a');
end

[TrialAligned, TrialInfo, ...
    ReplayAligned, ReplayInfo, ...
    TuningAligned, TuningInfo, ...
    AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

%MyUnits = 1:size(AllUnits.ChannelInfo,1); % all Units 
%MyUnits = [2 55 12 69 4 9 19 14 16 26 41 10]; %Q4
%MyUnits = [2 3 6 9 10 11 16 17 18 19 25 26 28 29 31 34 39 42 44 46 47 48 49 50 54 58 69 70 72 73 74 85 87 88 95 97]; % S12
MyUnits = [1 11 15 18 19 23 28 29 36 39 43 49 94]; %Q9

% get sniff time stamps and info for the sniffs we want to plot
[SelectedSniffs] = SelectSniffs_forKernelFits(TrialAligned, TrialInfo, [1 2 3]);

% % sort the sniffs
% sortorder = 1; % by Sniff Type, then sniff duration, then inh duration, then trial ID
% for whichodor = 1:3
%     SelectedSniffs{whichodor} = ...
%         SortSniffs(SelectedSniffs{whichodor}, sortorder);
% end
%%
PSTHbinsize = 10;

for unitcount = 1:numel(MyUnits)
    whichUnit = MyUnits(unitcount);
    
    SniffPSTHs = []; SniffParams = [];
    
    % get spike rasters
    for whichodor = 1:3
        % get the spike raster
        [SniffPSTHs_temp,SniffParams_temp] = ...
            SniffAlignedFRs(SelectedSniffs{whichodor}, AllUnits.Spikes{whichUnit}, ...
            'psthtype', 'binned', 'binsize', PSTHbinsize);
        
        SniffDim = size(SniffPSTHs_temp);
        SniffPSTHs(end+(1:SniffDim(1)),1:SniffDim(2)) = SniffPSTHs_temp;
        SniffParams = vertcat(SniffParams, SniffParams_temp);
    end
    % sniff params
    % 1-4: [currsniffstate currsniffloc currinhend currsniffend ...
    % 5-6:  currsniffTrialID currsniffIndex ...
    % 7:10: prevsniffstate prevsniffloc previnhstart previnhend]
    
    % starting kernels
    [StartingKernels] = InitialKernelEstimates(SniffPSTHs, SniffParams, ...
        'kernellength', 700, 'binsize', PSTHbinsize);
    
    %kernelsIn{unitcount} = StartingKernels;
    [kernelsout{unitcount},resnorm,residual,exitflag,output] = GetSniffKernels(StartingKernels, SniffParams, SniffPSTHs, ...
                 'binsize', PSTHbinsize);
    
%              fprintf(['took %d evals and has residual norm %f.\n'],...
%         output.funcCount, resnorm);
%     
%     [kernelsout2{unitcount},resnorm2,residual2,exitflag2,output2] = GetSniffKernels(kernelsout{unitcount}, SniffParams, SniffPSTHs, ...
%                  'binsize', PSTHbinsize);
%              
%     fprintf(['took %d evals and has residual norm %f.\n'],...
%         output2.funcCount, resnorm2);
    
end


%% compare PSTH
%MyUnits = sort(MyUnits);
for unitcount = 1:numel(MyUnits)
    whichUnit = MyUnits(unitcount);
    
    SniffPSTHs = []; SniffParams = [];
    
    % get spike rasters
    for whichodor = 1:3
        % get the spike raster
        [SniffPSTHs_temp,SniffParams_temp] = ...
            SniffAlignedFRs(SelectedSniffs{whichodor}, AllUnits.Spikes{whichUnit}, ...
            'psthtype', 'binned', 'binsize', PSTHbinsize);
        
        SniffDim = size(SniffPSTHs_temp);
        SniffPSTHs(end+(1:SniffDim(1)),1:SniffDim(2)) = SniffPSTHs_temp;
        SniffParams = vertcat(SniffParams, SniffParams_temp);
    end
    
    % predicted PSTHs
    [PSTHOut] = GetPredictedSniffPSTH(SniffParams,kernelsout{unitcount},'binsize', PSTHbinsize);
    
    %PSTHOut = PSTHOut*10;
    
    figure('Name',['unit ',num2str(whichUnit)],'Color','none');
    colormap(brewermap([100],'*YlGnBu'));
    odorstates = [-1 1 2 3]; 
    
    % plot stuff - sort sniffs by duration
    alims = [0 0]; blims = [0 0];
    for i = 1:4
        whichsniffs = find(SniffParams(:,1)==odorstates(i));
        
        % sort sniffs by duration
        [~,s] = sort(SniffParams(whichsniffs,4));
        
        subplot(2,8,(i*2)-1);
        imagesc(flipud(SniffPSTHs(whichsniffs(s),2:end)));
        alims = max([alims; get(gca,'CLim')]);
        
        subplot(2,8,(i*2)-0);
        imagesc(flipud(PSTHOut(whichsniffs(s),1:end)));
        blims = max([blims; get(gca,'CLim')]);
        
        % sort sniffs by odor location
        [~,s] = sort(SniffParams(whichsniffs,2));
        
        subplot(2,8,8 +(i*2)-1);
        imagesc(flipud(SniffPSTHs(whichsniffs(s),2:end)));
        
        subplot(2,8,8 + (i*2)-0);
        imagesc(flipud(PSTHOut(whichsniffs(s),1:end)));
    end
    for i = 2:2:16
        subplot(2,8,i-1);
        set(gca,'Clim',[0 alims(2)]);
        subplot(2,8,i);
        set(gca,'Clim',[0 blims(2)]);
    end
    set(gcf,'Position',[143 403 1653 420]);
    saveas(gcf,fullfile(FigPath,['unit ',num2str(whichUnit),'.png']));
end

%%
cf = []; 
for i = 1:numel(MyUnits) 
    cf(i,:) = kernelsout{i}([1 end]); 
    %cf(i,:) = kernelsout{i}([1 end-2:end]); 
end
figure('Name','location-coeffs & baselines')
subplot(2,1,1);
%plot(cf(:,2:4),'-o'); %,'XTicklabel',MyUnits);
plot(cf(:,2),'-o'); 
xticks(1:numel(MyUnits));
xticklabels(num2str(MyUnits'));

subplot(2,1,2);
plot(cf(:,1),'-o'); %,'XTicklabel',MyUnits);
xticks(1:numel(MyUnits));
xticklabels(num2str(MyUnits'));
saveas(gcf,fullfile(FigPath,['fittedcoeffs.png']));
