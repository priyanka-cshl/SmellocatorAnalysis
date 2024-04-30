% load the workspace for fitted kernels
load ('/home/priyanka/Desktop/sniffPSTHPredictions/S12/S12_20230731_sniffs.mat');

SelectedUnits = [29 25];

for unitcount = 1:numel(SelectedUnits)
    whichUnit = SelectedUnits(unitcount);
    
    unitID = find(MyUnits==whichUnit);
    
    % get the spike rasters
    SniffPSTHs = []; SniffParams = [];    
    for whichodor = 1:3
        [SniffPSTHs_temp,SniffParams_temp] = ...
            SniffAlignedFRs(SelectedSniffs{whichodor}, AllUnits.Spikes{whichUnit}, ...
            'psthtype', 'binned', 'binsize', PSTHbinsize);
        
        SniffDim = size(SniffPSTHs_temp);
        SniffPSTHs(end+(1:SniffDim(1)),1:SniffDim(2)) = SniffPSTHs_temp;
        SniffParams = vertcat(SniffParams, SniffParams_temp);
    end
    
    % get its predicted PSTH
    [PSTHOut] = GetPredictedSniffPSTH(SniffParams,kernelsout{unitID},'binsize', PSTHbinsize);
    
    % plotting
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
    
    dt = 0.010; % in s
    PSTHtemp = PSTHOut(whichsniffs(s),1:end);
    figure; hold on
    for s = 1:size(PSTHtemp,1)
        [T] = PechePourPoisson(100*PSTHtemp(s,:),dt);
        plot(T,s + 0*T, '.k','Markersize', 0.5);
    end
    
end

% 
%     %%
%     
%     SniffPSTHs = []; SniffParams = [];
%     
%     % get spike rasters
%     for whichodor = 1:3
%         % get the spike raster
%         [SniffPSTHs_temp,SniffParams_temp] = ...
%             SniffAlignedFRs(SelectedSniffs{whichodor}, AllUnits.Spikes{whichUnit}, ...
%             'psthtype', 'binned', 'binsize', PSTHbinsize);
%         
%         SniffDim = size(SniffPSTHs_temp);
%         SniffPSTHs(end+(1:SniffDim(1)),1:SniffDim(2)) = SniffPSTHs_temp;
%         SniffParams = vertcat(SniffParams, SniffParams_temp);
%     end
%     
%     % predicted PSTHs
%     [PSTHOut] = GetPredictedSniffPSTH(SniffParams,kernelsout{unitcount},'binsize', PSTHbinsize);
%     
%     %PSTHOut = PSTHOut*10;
%     
%     
%     set(gcf,'Position',[143 403 1653 420]);
%     saveas(gcf,fullfile(FigPath,['unit ',num2str(whichUnit),'.png']));
% end