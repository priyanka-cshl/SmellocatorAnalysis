% function to plot spikes aligned by trial types
function [] = PlotLocationOffsetTraces(Traces,TrialInfo)

global MyFileName
% global subplotcol
% if subplotcol == 0
%     subplotcol = 5;
startoffset = 1;

% end
global SampleRate;
plotX = 1; % one for each odor and one pooled
plotY = 2; % behavior, aligned to perturbation start (raster, psth), behavior, aligned to feedback start (raster, psth)

xbins = -2000:1:2000;

X = colormap(brewermap([11],'*RdBu'));

% [Traces, TrialInfo, TargetZones] = ParseTrials(MyData, MySettings, TargetZones, sessionstart, sessionstop);
% [spiketimes] = Spikes2Trials(myephysdir);

if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    u = unique(TrialInfo.Perturbation(x));
    perturbation_params = cell2mat(TrialInfo.Perturbation(x,2));
end


% split by different odors
for whichodor = 4
    
    if whichodor<4
        % extract trials that were perturbed
        %whichtrials = find(TrialInfo.Odor==whichodor & (TrialInfo.Perturbation(:,1)==6 | TrialInfo.Perturbation(:,1)==7));
        whichtrials = intersect(x,find(TrialInfo.Odor==whichodor));
    else
        %whichtrials = find(TrialInfo.Perturbation(:,1)==6 | TrialInfo.Perturbation(:,1)==7);
        whichtrials = x;
    end
    
    % Re-Sort to differenciate offset sizes
    if ~isempty(whichtrials)
        order = [find(perturbation_params(:,1)>121); find(perturbation_params(:,1)<121)];
        whichtrials = whichtrials(order);
        perturbation_params = perturbation_params(order,:);
        div_line = find(perturbation_params(:,1)<121,1,'first');
    end
    
    if ~isempty(whichtrials)
        
        % each trial
        for trial = 1:numel(whichtrials)
            trial_idx = whichtrials(trial); % trial ID
            perturbation_idx = trial;
            lever = cell2mat(Traces.Lever(trial_idx)); % in samples @500 Hz
            
            perturbationstart = startoffset+perturbation_params(perturbation_idx,2); % in seconds w.r.t. trace start
            feedbackstart = startoffset+perturbation_params(perturbation_idx,3); % in seconds w.r.t. trace start
            
            thisTrialRewards = cell2mat(TrialInfo.Reward(trial_idx));
            if ~isempty(thisTrialRewards) && any(thisTrialRewards>0)
                reward = thisTrialRewards(find(thisTrialRewards>0,1,'first')); % in seconds w.r.t. trace start
            else
                reward = NaN;
            end
            
            %% 1: align to perturbation start
            
            % plot the lever traces
            subplot(1,2,1);
            hold on
            time_idx = (1000/SampleRate)*(1:numel(lever)) - round(perturbationstart*1000);
            if trial < div_line
                plot(time_idx,lever,'color',X(2,:));
            else
                plot(time_idx,lever,'color',X(end-1,:));
            end
            
            % overlay feedback start (as red ticks)
            feedbacktick = round((feedbackstart - perturbationstart)*1000);
            line([feedbacktick feedbacktick],...
                [trial-1 trial],'Color','r', 'LineWidth', 2);
            
            %% 2: align to feedback start

            % plot the lever traces
            subplot(1,2,2);
            hold on
            time_idx = (1000/SampleRate)*(1:numel(lever)) - round(feedbackstart*1000);
            if trial < div_line
                plot(time_idx,lever,'color',X(2,:));
            else
                plot(time_idx,lever,'color',X(end-1,:));
            end
            
            % overlay feedback start (as red ticks)
            if ~isnan(reward)
                rewardtick = round((reward - feedbackstart)*1000);
                line([rewardtick rewardtick],...
                    [trial-1 trial],'Color','g', 'LineWidth', 2);
            end
           
        end

        subplot(1,2,1);
        line([0 0],[-0.5 5.5],'Color','b', 'LineWidth', 1);
        set(gca,'XLim',[xbins(1) xbins(end)],'YLim',[-0.5 5.5]);
        set(gca,'YTick',[],'XTick',[xbins(1) 0 xbins(end)]);
        
        subplot(1,2,2);
        line([0 0],[-0.5 5.5],'Color','r', 'LineWidth', 1);
        set(gca,'XLim',[xbins(1) xbins(end)],'YLim',[-0.5 5.5]);
        set(gca,'YTick',[],'XTick',[xbins(1) 0 xbins(end)]);
        
    end
end

end

