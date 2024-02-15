function [AllOdorSniffs] = SelectSniffs(TrialAligned, TrialInfo, OdorList, varargin)

%%narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('includeITI', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('includePerturbations', false, @(x) islogical(x) || x==0 || x==1);
% params.addParameter('BehaviorSampleRate', 500, @(x) isnumeric(x)); % 0 - sniff duration, 1 - inhalation duration
% params.addParameter('sortorder', 0, @(x) isnumeric(x)); % 0 - sniff duration, 1 - inhalation duration
% params.addParameter('whichunits', [], @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
includeITI = params.Results.includeITI;

AllOdorSniffs = [];

for n = 1:numel(OdorList)
    
    whichodor = OdorList(n);
    AllSniffs = [];
    
    %% for closed loop or replay trials
    % get all unpertubed trials
    whichtrials = find(cellfun(@isempty, TrialInfo.Perturbation(:,1))); % vanilla closed-loop
    % also include trials marked as OL-Template - they are also just control trials
    whichtrials = sort([whichtrials; find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'))]); % also used the template trials (no-perturbations)
    % only keep trials of a particular odor
    if isfield(TrialInfo, 'Odor')
        whichtrials(find(TrialInfo.Odor(whichtrials)~=whichodor),:) = [];
    end
    
    for i = 1:numel(whichtrials)
        thisTrial       = whichtrials(i); % Trial ID
        thisTrialStart  = TrialAligned.RefTS(thisTrial,1);
        
        if ~includeITI
            whichsniffs = find(TrialAligned.Sniffs{thisTrial}(:,2)==odorstate);
        else
            whichsniffs = (1:size(TrialAligned.Sniffs{thisTrial},1))';
        end
        SniffTimestamps = TrialAligned.Sniffs{thisTrial}(whichsniffs,:);
        SniffTimestamps(:,5:11) = SniffTimestamps(:,5:11) + thisTrialStart;
        AllSniffs       = [AllSniffs; SniffTimestamps];
    end
    
    % add a column for current sniff and inhalation duration
    AllSniffs(:,15:16) = AllSniffs(:,[9 8]) - AllSniffs(:,7);
    
    AllOdorSniffs{whichodor} = AllSniffs;
    
    % AllSniffs = sortrows(AllSniffs, [2 4 15]); % odorstate, snifftype and duration
end
end

% % get halt trials if any
% % also collect perturbation trials
% perturbationTrials = intersect(find(strncmpi(TrialInfo.Perturbation(:,1),'Halt-Flip',9)), ...
%     find(TrialInfo.Odor==whichodor));
% 
% if numel(perturbationTrials) < 2
%     perturbationTrials = [];
% end
% 
% % assemble a list of sniffs 
% % [inh-start inh-end next-inh sniffID snifftype sniff-duration inh-duration Trial ID]
% 
% AllSniffs = [];
% for s = 1:numel(allTrials) % every trial
%     thisTrialSniffs = TrialAlignedSniffs{allTrials(s)}; 
%     thisTrialSniffs(:,15) = thisTrialSniffs(:,3) - thisTrialSniffs(:,1); % col 15 = sniff duration
%     thisTrialSniffs(:,16) = thisTrialSniffs(:,2) - thisTrialSniffs(:,1); % col 16 = inhalation duration
%     thisTrialSniffs(:,17) = allTrials(s); % col 17 = trial ID
%     AllSniffs = vertcat(AllSniffs, thisTrialSniffs);
% end
% 
% % add perturbation sniffs if any
% if ~isempty(perturbationTrials)
%     for s = 1:numel(perturbationTrials)
%         thisTrialSniffs = TrialAlignedSniffs{perturbationTrials(s)}; 
%         haltperiod = TrialInfo.Perturbation{perturbationTrials(s),2}(1:2)./BehaviorSampRate;
%         haltsniffs = intersect(find(thisTrialSniffs(:,1)>=haltperiod(1)),find(thisTrialSniffs(:,1)<=haltperiod(2)));
%         thisTrialSniffs = thisTrialSniffs(haltsniffs,:);
%         thisTrialSniffs(:,15) = thisTrialSniffs(:,3) - thisTrialSniffs(:,1);
%         thisTrialSniffs(:,16) = thisTrialSniffs(:,2) - thisTrialSniffs(:,1);
%         thisTrialSniffs(:,17) = perturbationTrials(s);
%         % change snifftype 
%         thisTrialSniffs(:,5) = 4;
%         AllSniffs = vertcat(AllSniffs, thisTrialSniffs);
%     end
%     
%     haltlocation = TrialInfo.Perturbation{perturbationTrials(2),2}(3);
%     
%     % also pull out control sniffs that share the same location
%     ctrl_sniffs = intersect(find((AllSniffs(:,5)>-1)&(AllSniffs(:,5)<4)), ...
%                     find(round(AllSniffs(:,6)/10)==haltlocation/10));
%     AllSniffs(ctrl_sniffs,5) = 3;
%     
% else
%     haltlocation = NaN;            
% end
% 
% if haltmode
%     AllSniffs(find(AllSniffs(:,5)<3),:) = []; % only keep closeloop sniffs that share the halt location
% end
% 
% switch sortby
%     case 0
%         % sort sniff List by Sniff Type, then sniff duration, then inh duration, then trial ID
%         AllSniffs = sortrows(AllSniffs,[5 15 16 17]);
%     case 1
%         % sort sniff List by Sniff Type, then inh duration, then sniff duration, then trial ID
%         AllSniffs = sortrows(AllSniffs,[5 16 15 17]);
%     case 2
%         AllSniffs(find(AllSniffs(:,5)==0),5) = 1; 
%         AllSniffs(find(AllSniffs(:,5)==1),5) = 1;
%         AllSniffs(find(AllSniffs(:,5)==2),5) = 1; 
%         % sort by sniff type, then odor location (col 6), then sniff duration, then inh duration, then trial ID
%         AllSniffs = sortrows(AllSniffs,[5 6 16 15 17]);
% end