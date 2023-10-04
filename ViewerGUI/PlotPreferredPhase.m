function [] = PlotPreferredPhase(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo)

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
whichUnit = abs(whichUnit);

% hack to prevent OL-Template trials to be considered as perturbed trials
f = find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'));
if ~isempty(f)
    for i = 1:numel(f)
        TrialInfo.Perturbation{f(i),1} = [];
    end
end

thisUnitSpikes = AlignedSpikes(:,whichUnit);
whichodor = whichOdor;
% get the trials for a given odor
whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
    find(TrialInfo.Odor==whichodor));
allTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) TrialInfo.Duration(whichTrials)]; %#ok<AGROW>

if plotting
    B = []; T = [];
    for x = 1:size(whichTrials,1) % every trial
        thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1}; % all spikes in this Trial
        
        % baseline sniffs - all spikes less than 0
        trialSpikes = thisTrialSpikes(thisTrialSpikes<floor(Events(whichTrials(x,1),2)));
        baselineSpikes = trialSpikes(trialSpikes<0);
        trialSpikes(trialSpikes<1) = [];
        
        B = horzcat(B,baselineSpikes - floor(baselineSpikes));
        T = horzcat(T,trialSpikes - floor(trialSpikes));
    end
    
    histogram(B,[0:0.05:1],'Normalization','probability');
    hold on
    histogram(T,[0:0.05:1],'Normalization','probability');
end

end