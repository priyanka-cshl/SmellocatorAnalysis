%% load GlomSession into work space
clear

load('C:\Users\blomk\OneDrive - UvA\Internship 2\Data\bulb_imaging\testmouse2_black\12-Oct-2022_1\AllGloms.mat');

% get stimulus period
preStim = GlomSession.TrialSequence(1,5);
stimulus = (GlomSession.TrialSequence(1,7) - GlomSession.TrialSequence(1,6) + 1);
endStim = preStim + stimulus;
postStim = 160;

% remove dark frames
GlomSession.Traces = GlomSession.Traces(:,1:(endStim + postStim),:);

% baseline correct traces with each glom as its own baseline
for trial = 1:length(GlomSession.TrialSequence)
    
    for ROI = 1:length(GlomSession.ROI_index)
        
        baseline = mean(GlomSession.Traces(ROI,1:(GlomSession.TrialSequence(1,5)), trial));
        
        % dF/F
        GlomSession.Traces(ROI, :, trial) = GlomSession.Traces(ROI, :, trial) - baseline;
        GlomSession.Traces(ROI, :, trial) = GlomSession.Traces(ROI, :, trial)/baseline;
    end
    
end

% pool all glomeruli
avg_allGlom = squeeze(mean(GlomSession.Traces, 1));

%% responsive glomeruli

% GlomSession.Traces previously baseline corrected
% indentify which glomerulus was actived by which odor

% baseline correct traces with each glom as its own baseline
[trials] = find((GlomSession.TrialSequence(:,2) == 1));
%(GlomSession.TrialSequence(:,1) == 0) & 
responsive_ROIs = zeros(length(trials), (length(GlomSession.ROI_index)));

for trial = 1:length(trials)
    
    for ROI = 1:length(GlomSession.ROI_index)
        
        % all baseline frames
        baseline = GlomSession.Traces(ROI,(1: preStim), trial);
        % all stimulus frames
        odorOn = GlomSession.Traces(ROI,((preStim + 1) : 160), trial);
        
        % wilcoxon signed rank test if odor + 2SD > baseline
        [p,h] = signrank(baseline, odorOn, 'alpha', 0.001, 'tail', 'left');
        %(odorOn + (10 * std(odorOn)))
        
        if h == 0
            responsive_ROIs(trial, GlomSession.ROI_index(ROI)) = h;
            plot(GlomSession.Traces(ROI,:, trial))
            pause;
        end
        
    end
end


% [nr_responses, rois] = find(responsive_ROIs == 1);
% air_resp = length(unique(rois));
% remove any glom responsive to air
 
%% combine data from 2 animals

TuningCurves_2 = load('C:\Users\blomk\OneDrive - UvA\Internship 2\Data\bulb_imaging\TuningCurves_2.mat');
TuningCurves_2 = TuningCurves_2.TuningCurves;
TuningCurves_3 = load('C:\Users\blomk\OneDrive - UvA\Internship 2\Data\bulb_imaging\TuningCurves_3.mat');
TuningCurves_3 = TuningCurves_3.TuningCurves;

TuningCurves = cat(1, TuningCurves_2, TuningCurves_3);

% get stimulus period
preStim = 80;
stimulus = 20; % high: 80
endStim = preStim + stimulus;
postStim = 160;


%% PSTH per odor
% average normalized glom response
%avg_allGlom = squeeze(mean(GlomSession.Traces, 1));

figure
% plot response per odor
for j = 1:length(GlomSession.Odors)
    
    % find only center location
    [trials] = find((GlomSession.TrialSequence(:,2) == j) & (GlomSession.TrialSequence(:,1) == 0));
    
    % plot mean figures
    %figure
    %hold on
    subplot(2, 2, j)
    
    % std
    shadedErrorBar([],(mean(avg_allGlom(:, trials), 2)), (std(avg_allGlom(:,trials)')), 'b');
    
    % SEM
    % shadedErrorBar([],(mean(avg_allGlom(:, trials), 2)),(std(avg_allGlom(:,trials)') / sqrt(length(trials))), 'r')
    
    title(sprintf('Odor %i, location 0', j));
    xlabel('frame nr');
    xline(preStim)
    xline(endStim)
    %ylabel('raw data');
    %legend('mean trace');
    
end

%% each glom's tuning curve, with average
%TuningCurves = zeros(size(GlomSession.Traces,1),length(GlomSession.Locations),length(GlomSession.Odors));

figure;

% plot for each odor
for m = 1:length(GlomSession.Odors)
    
    subplot(2, 2, m)
    % as a function of location
    for n = 1:length(GlomSession.Locations)
        
        % find loc for specific odor
        [trials] = find((GlomSession.TrialSequence(:,2) == m) & (GlomSession.TrialSequence(:,1) == GlomSession.Locations(n)));
        
        for p = 1:size(GlomSession.Traces,1)
            TuningCurves(p,n,m) = mean(mean(GlomSession.Traces(p, preStim:endStim, trials)));
        end
    end
    
    % single tuning curves
    plot(GlomSession.Locations, TuningCurves(:, :, m), 'k:.');
    hold on
    % average tuning curve
    plot(GlomSession.Locations, (mean(TuningCurves(:, :, m))), 'c.-', 'LineWidth', 1.25);
    
    % ylim([-0.1 0.9])
    xticks(GlomSession.Locations)
    title(sprintf('Odor %i', m));
    xlabel('odor space');
    
    
end

%save('TuningCurves_3.mat', 'TuningCurves');

%% plot response of a few glom for all 4 odors

for i = 1:length(GlomSession.Odors)
    subplot(2,2,i)
    hold on
    plot(GlomSession.Locations, TuningCurves(2,:,i));
    plot(GlomSession.Locations, TuningCurves(12,:,i));
    plot(GlomSession.Locations, TuningCurves(20,:,i));
    plot(GlomSession.Locations, TuningCurves(7,:,i)); %45
    plot(GlomSession.Locations, TuningCurves(15,:,i)); % 93
    plot(GlomSession.Locations, TuningCurves(74,:,i));
    title(sprintf('Odor %i', i));
    %xticks(GlomSession.Locations([1 4 7 10 13]))
end


%% population activity during odor period
figure
% plot for each odor
for m = 1:length(GlomSession.Odors)
    
    %figure
    subplot(2, 2, m)
    
    % as a function of location
    for n = 1:length(GlomSession.Locations)
        
        % find loc for specific odor
        [trials] = find((GlomSession.TrialSequence(:,2) == m) & (GlomSession.TrialSequence(:,1) == GlomSession.Locations(n)));
        
        hold on
        
        % replot this so that it will plot in 1 go, see section below
        plot(GlomSession.Locations(n), (mean(mean(avg_allGlom(preStim:endStim, trials)))), '.');
        xticks(GlomSession.Locations)
        ylim([-0.1 0.25])
        title(sprintf('Odor %i', m));
        xlabel('odor space');
        
    end
end



%% to remember
%GlomSession.Traces(1,:,1); % plot 1 trace
%avg_allTrials = mean(GlomSession.Traces, 3);
%plot(avg_allTrials') % 1 trace = 1 glom
%avg_all = mean(avg_allTrials);

%Frames = length(GlomSession.Traces(1,:,1));
%baseline_frames = GlomSession.TrialSequence(1,5);

% Baseline = mean(MyStack(:,:,:,BaselineFrames),4);
% Stimulus = mean(MyStack(:,:,:,StimulusFrames),4);
% Ratio(:,:,thisTrial) = (Stimulus - Baseline)./Baseline;

%r = (find(GlomSession.TrialSequence(:,2) == 1) && find(GlomSession.TrialSequence(:,1) == 0))
% meanTrace = mean(avg_allGlom(:,trialNr), 2);
%SEM = std(avg_allGlom(:,trialNr)')/sqrt(length(trialNr));

%plot((mean(avg_allGlom(:, trialNr), 2)), 'k');
%plot(mean(avg_allGlom(:, trialNr), 2)' + (std(avg_allGlom(:,trialNr)') / sqrt(length(trialNr))));
%plot(mean(avg_allGlom(:, trialNr), 2)' - (std(avg_allGlom(:,trialNr)') / sqrt(length(trialNr))));

%errorbar((mean(mean(avg_allGlom(preStim:endStim, trials)))), GlomSession.Locations(n), (mean(std(avg_allGlom(preStim:endStim, trials)))));

% %% deduct average baseline intensity from entire trace
% for i = 1:length(GlomSession.TrialSequence)
%
%     avg_allGlom(:, i) = avg_allGlom(:, i) - mean(avg_allGlom(1:((GlomSession.TrialSequence(1,5)-10)), i));
%
%     % Ratio(:,:,thisTrial) = (Stimulus - Baseline)./Baseline;
% end
