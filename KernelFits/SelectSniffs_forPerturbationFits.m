function [AllPerturbationSniffs] = SelectSniffs_forPerturbationFits(TrialAligned, TrialInfo) %, OdorList)

AllPerturbationSniffs = [];

% for n = 1:numel(OdorList)
%     
%     whichodor = OdorList(n);
    AllSniffs = [];
    
    %% for closed loop or replay trials
    if ~isfield(TrialInfo,'TuningTrialID')
        % get all unpertubed trials
        whichtrials = find(cellfun(@isempty, TrialInfo.Perturbation(:,1))); % vanilla closed-loop
        % also include trials marked as OL-Template - they are also just control trials
        whichtrials = sort([whichtrials; find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'))]); % also used the template trials (no-perturbations)
    else
        whichtrials = (1:size(TrialInfo.Odor,1))';
    end
    % perturbed trials
    AllTrials = 1:size(TrialInfo.TrialID,1);
    perturbed = AllTrials(~ismember(AllTrials,whichtrials));
    whichtrials = perturbed;
    
%     % only keep trials of a particular odor
%     if isfield(TrialInfo, 'Odor')
%         whichtrials(find(TrialInfo.Odor(whichtrials)~=whichodor),:) = [];
%     end
    
    for i = 1:numel(whichtrials)
        thisTrial       = whichtrials(i); % Trial ID
        %whichsniffs = (2:size(TrialAligned.Sniffs{thisTrial},1))'; % ignore the first sniff after trial OFF - uncertainty about the prev sniff state
        
        if regexp(TrialInfo.Perturbation{thisTrial,1},'Halt')
            perturbation_duration = TrialInfo.Perturbation{thisTrial,2}(2) - TrialInfo.Perturbation{thisTrial,2}(1); % in indices
            perturbation_duration = perturbation_duration/500; % in seconds for behavior sampling rate of 500hz
            t1 = TrialAligned.Events(thisTrial,5); % halt start
            t2 = t1 + perturbation_duration; %halt end
            % only keep sniffs for whom inhaltion started after
            % perturbation start and inhalation ended before perturbation
            % end
            whichsniffs = intersect(find(TrialAligned.Sniffs{thisTrial}(:,6)>=t1), ...
                find(TrialAligned.Sniffs{thisTrial}(:,7)<=t2) );
            
        else
            whichsniffs = (2:size(TrialAligned.Sniffs{thisTrial},1))'; % ignore the first sniff after trial OFF - uncertainty about the prev sniff state
        end
        thisTrialStart  = TrialAligned.RefTS(thisTrial,1); % actual trial start time in behavior timebase
        SniffTimestamps = TrialAligned.Sniffs{thisTrial}(whichsniffs,:);
        SniffTimestamps(:,5:11) = SniffTimestamps(:,5:11) + thisTrialStart;
        SniffTimestamps(:,18)   = TrialAligned.Sniffs{thisTrial}(whichsniffs-1,2); % previous sniff's Odor/Air State
        AllSniffs       = [AllSniffs; SniffTimestamps];
        
    end
    
    % add a column for current sniff and inhalation duration
    AllSniffs(:,15:16) = AllSniffs(:,[9 8]) - AllSniffs(:,7);
    % add a column for previous sniff duration
    AllSniffs(:,17) = AllSniffs(:,7) - AllSniffs(:,5);
    
    
    AllPerturbationSniffs = AllSniffs;
    
    % AllSniffs = sortrows(AllSniffs, [2 4 15]); % odorstate, snifftype and duration
% end
end