% function to plot spikes aligned by trial types
function [] = PlotLocationOffsetPSTH_2018(whichcluster,Traces,TrialInfo,spiketimes)

global MyFileName
% global subplotcol
% if subplotcol == 0
%     subplotcol = 5;
% end
global SampleRate;
plotX = 4; % one for each odor and one pooled
plotY = 6; % behavior, aligned to perturbation start (raster, psth), behavior, aligned to feedback start (raster, psth)

xbins = -2000:1:2000;

X = colormap(brewermap([11],'*RdBu'));

% [Traces, TrialInfo, TargetZones] = ParseTrials(MyData, MySettings, TargetZones, sessionstart, sessionstop);
% [spiketimes] = Spikes2Trials(myephysdir);

if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    u = unique(TrialInfo.Perturbation(x));
end


% split by different odors
for whichodor = 1:4
    
    if whichodor<4
        % extract trials that were perturbed
        %whichtrials = find(TrialInfo.Odor==whichodor & (TrialInfo.Perturbation(:,1)==6 | TrialInfo.Perturbation(:,1)==7));
        whichtrials = intersect(x,find(TrialInfo.Odor==whichodor));
    else
        %whichtrials = find(TrialInfo.Perturbation(:,1)==6 | TrialInfo.Perturbation(:,1)==7);
        whichtrials = x;
    end
    
    % Re-Sort to differenciate offset sizes
    whichtrials = [whichtrials(TrialInfo.Perturbation(whichtrials,2)>121); whichtrials(TrialInfo.Perturbation(whichtrials,2)<121)];
    div_line = find(TrialInfo.Perturbation(whichtrials,2)<121,1,'first');
    
    if ~isempty(whichtrials)
        psth = zeros(4,numel(xbins)); % aligned to perturbation start
        psth2 = zeros(4,numel(xbins)); % aligned to feedback start
        
        % each trial
        for trial = 1:numel(whichtrials)
            trial_idx = whichtrials(trial); % trial ID
            thisTrialSpikes = cell2mat(spiketimes(whichcluster).spikes(TrialInfo.TrialID(trial_idx)));
            lever = cell2mat(Traces.Lever(trial_idx)); % in samples @500 Hz
            
            perturbationstart = TrialInfo.PerturbationStart(trial_idx); % in seconds w.r.t. trace start
            feedbackstart = TrialInfo.FeedbackStart(trial_idx); % in seconds w.r.t. trace start
            
            thisTrialRewards = cell2mat(TrialInfo.Reward(trial_idx));
            if ~isempty(thisTrialRewards) && any(thisTrialRewards>0)
                reward = thisTrialRewards(find(thisTrialRewards>0,1,'first')); % in seconds w.r.t. trace start
            else
                reward = NaN;
            end
            
            %% 1: align to perturbation start
            
            % plot the lever traces
            subplot(plotX,plotY,whichodor*plotY - plotY + 1);
            hold on
            time_idx = (1000/SampleRate)*(1:numel(lever)) - round(perturbationstart*1000);
            if trial < div_line
                plot(time_idx,lever,'color',X(2,:));
            else
                plot(time_idx,lever,'color',X(end-1,:));
            end
            
            % plot spikes
            subplot(plotX,plotY,whichodor*plotY - plotY + 2);
            hold on
            myspikes = round((thisTrialSpikes - perturbationstart)*1000); % in ms
            for eachspike = 1:numel(myspikes) % plot raster line
                line([myspikes(eachspike) myspikes(eachspike)],...
                    [trial-1 trial],'Color','k');
                hold on
            end
            
            % overlay feedback start (as red ticks)
            feedbacktick = round((feedbackstart - perturbationstart)*1000);
            line([feedbacktick feedbacktick],...
                [trial-1 trial],'Color','r', 'LineWidth', 2);
            
            % update psth
            for inx = 1:length(xbins)-1 %calculate # spikes in ms bins
                psth(2-(trial<div_line),inx) = psth(2-(trial<div_line),inx) + ...
                    numel(find(myspikes>xbins(inx) & myspikes<=xbins(inx+1)));
            end
            lowlim = max(-round(perturbationstart*1000),xbins(1));
            uplim = min(round((TrialInfo.Timestamps(trial_idx,2)-perturbationstart)*1000),xbins(end));
            psth(4-(trial<div_line),find(xbins==lowlim):find(xbins==uplim)) = psth(4-(trial<div_line),find(xbins==lowlim):find(xbins==uplim)) + 1;
            
            
            %% 2: align to feedback start

            % plot the lever traces
            subplot(plotX,plotY,whichodor*plotY - plotY + 4);
            hold on
            time_idx = (1000/SampleRate)*(1:numel(lever)) - round(feedbackstart*1000);
            if trial < div_line
                plot(time_idx,lever,'color',X(2,:));
            else
                plot(time_idx,lever,'color',X(end-1,:));
            end
            
            % plot spikes
            subplot(plotX,plotY,whichodor*plotY - plotY + 5);
            hold on
            myspikes = round((thisTrialSpikes - feedbackstart)*1000); % in ms
            for eachspike = 1:numel(myspikes) % plot raster line
                line([myspikes(eachspike) myspikes(eachspike)],...
                    [trial-1 trial],'Color','k');
                hold on
            end
            
            % overlay feedback start (as red ticks)
            if ~isnan(reward)
                rewardtick = round((reward - feedbackstart)*1000);
                line([rewardtick rewardtick],...
                    [trial-1 trial],'Color','g', 'LineWidth', 2);
            end
            % update psth
            for inx = 1:length(xbins)-1 %calculate # spikes in ms bins
                psth2(2-(trial<div_line),inx) = psth2(2-(trial<div_line),inx) + ...
                    numel(find(myspikes>xbins(inx) & myspikes<=xbins(inx+1)));
            end
            lowlim = max(-round(feedbackstart*1000),xbins(1));
            uplim = min(round((TrialInfo.Timestamps(trial_idx,2)-feedbackstart)*1000),xbins(end));
            psth2(4-(trial<div_line),find(xbins==lowlim):find(xbins==uplim)) = psth2(4-(trial<div_line),find(xbins==lowlim):find(xbins==uplim)) + 1;
            
        end
        
        % Normalize psths
        for bin = 1:numel(xbins)
            psth(1,bin) = psth(1,bin)/psth(3,bin); % divide by trial counts
            psth(2,bin) = psth(2,bin)/psth(4,bin); % divide by trial counts
            psth2(1,bin) = psth2(1,bin)/psth2(3,bin); % divide by trial counts
            psth2(2,bin) = psth2(2,bin)/psth2(4,bin); % divide by trial counts
        end
        psth(1,:) = 1000*psth(1,:); % Hz
        psth(2,:) = 1000*psth(2,:); % Hz
        sigma = 20;
        conv_psth(1,:) = convPSTH(psth(1,:), sigma);
        conv_psth(2,:) = convPSTH(psth(2,:), sigma);
        
        psth2(1,:) = 1000*psth2(1,:); % Hz
        psth2(2,:) = 1000*psth2(2,:); % Hz
        conv_psth2(1,:) = convPSTH(psth2(1,:), sigma);
        conv_psth2(2,:) = convPSTH(psth2(2,:), sigma);
        
        subplot(plotX,plotY,whichodor*plotY - plotY + 3);
        hold on
        plot(xbins, conv_psth(1,:), 'color', X(2,:), 'LineWidth', 1);
        plot(xbins, conv_psth(2,:), 'color', X(end-1,:), 'LineWidth', 1);
        line([0 0],get(gca,'YLim'),'Color','b', 'LineWidth', 1);
        set(gca,'XLim',[xbins(1) xbins(end)],'XTick',[xbins(1) 0 xbins(end)]);
        
        subplot(plotX,plotY,whichodor*plotY - plotY + 6);
        hold on
        plot(xbins, conv_psth2(1,:), 'color', X(2,:), 'LineWidth', 1);
        plot(xbins, conv_psth2(2,:), 'color', X(end-1,:), 'LineWidth', 1);
        line([0 0],get(gca,'YLim'),'Color','r', 'LineWidth', 1);
        set(gca,'XLim',[xbins(1) xbins(end)],'XTick',[xbins(1) 0 xbins(end)]);
        
        
        subplot(plotX,plotY,whichodor*plotY - plotY + 2);
        line([0 0],[-1 trial+1],'Color','b', 'LineWidth', 1);
        if div_line
            line([xbins(1) xbins(end)],[div_line div_line],'Color','k','LineStyle',':', 'LineWidth', 1);
        end
        set(gca,'XLim',[xbins(1) xbins(end)],'YLim',[-1 trial+1]);
        set(gca,'YTick',[],'XTick',[xbins(1) 0 xbins(end)]);
        
        subplot(plotX,plotY,whichodor*plotY - plotY + 5);
        line([0 0],[-1 trial+1],'Color','b', 'LineWidth', 1);
        if div_line
            line([xbins(1) xbins(end)],[div_line div_line],'Color','k','LineStyle',':', 'LineWidth', 1);
        end
        set(gca,'XLim',[xbins(1) xbins(end)],'YLim',[-1 trial+1]);
        set(gca,'YTick',[],'XTick',[xbins(1) 0 xbins(end)]);
        
        subplot(plotX,plotY,whichodor*plotY - plotY + 1);
        line([0 0],[0 5],'Color','b', 'LineWidth', 1);
        set(gca,'XLim',[xbins(1) xbins(end)],'YLim',[0 5]);
        set(gca,'YTick',[],'XTick',[xbins(1) 0 xbins(end)]);
        
        subplot(plotX,plotY,whichodor*plotY - plotY + 4);
        line([0 0],[0 5],'Color','r', 'LineWidth', 1);
        set(gca,'XLim',[xbins(1) xbins(end)],'YLim',[0 5]);
        set(gca,'YTick',[],'XTick',[xbins(1) 0 xbins(end)]);
        
    end
end

end

