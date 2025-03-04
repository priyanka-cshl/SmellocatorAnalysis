%% mini function to recalculate sniff TS in behavior time base - trial by trial

function [SniffTimeStamps] = TrialWiseSniffs(TrialInfo,Traces,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('SDfactor', 2.5, @(x) isnumeric(x));
params.addParameter('dlgoverride', logical(0), @(x) islogical(x));

% extract values from the inputParser
params.parse(varargin{:});
SDfactor = params.Results.SDfactor;
nodlg = params.Results.dlgoverride;

%% sniff peak-valley detection
SniffTimeStamps = [];
for trialID = 1:numel(TrialInfo.Odor)
    % detect sniff timestamps per trial
    TraceTimeStamps = Traces.Timestamps{trialID};
    Thermistor      = Traces.Sniffs{trialID};
    OdorLocation    = Traces.OdorLocation{trialID};
    thisTrialSniffs = ProcessThermistorData([TraceTimeStamps Thermistor OdorLocation],'SDfactor',SDfactor,'dlgoverride',nodlg);

    % apppend to AllSniffs and also handle end cases
    if trialID > 1
        % if first sniff started earlier than this trace
        if isnan(thisTrialSniffs(1,1))
            [~,f] = min(sum(abs(SniffTimeStamps(:,2:3) - thisTrialSniffs(1,2:3)),2));
            if ~isempty(f)
                thisTrialSniffs(1,:) = SniffTimeStamps(f,1:5);
            end
        end
        % if the last sniff in Allsniffs was truncated
        if isnan(SniffTimeStamps(end,3))
            [~,f] = min(sum(abs(thisTrialSniffs(:,1:2) - SniffTimeStamps(end,1:2)),2));
            if ~isempty(f)
                SniffTimeStamps(end,1:5) = thisTrialSniffs(f,1:5);
            end
        end
        % flag overlapping sniffs
        %f = find(ismember(SniffTimeStamps(:,1:4),thisTrialSniffs(1,1:4),'rows'));
        f = find(SniffTimeStamps(:,1)>=thisTrialSniffs(1,1)-0.005);
        if ~isempty(f)
            SniffTimeStamps(f:end,7) = -SniffTimeStamps(f:end,7);
        end
    end
    
    % add trialstate? -1 - air off, 0,1,2,3 - air, odor1, odor2, odor3
    if isfield(TrialInfo,'SessionTimestamps')
        thisTrialStartStop = TrialInfo.SessionTimestamps(trialID,1:2);
        thisTrialOdorStart = TrialInfo.OdorStart(trialID,2); % in seconds before trial start
    else
        thisTrialStartStop = [find(diff(Traces.TrialState{trialID})>0,1,'first') find(diff(Traces.TrialState{trialID})<0,1,'last')];
        thisTrialStartStop = Traces.Timestamps{trialID}(thisTrialStartStop);
        thisTrialOdorStart = TrialInfo.OdorStart(trialID,1); % in seconds before trial start
    end
    
    thisTrialOdorStart = thisTrialStartStop(1) + thisTrialOdorStart; 
    thisTrialOdorStop  = thisTrialStartStop(2);
    
    thisTrialSniffs(:,6) = -1; % mark all as being trialOFF
    % pull up all inhalations that happened after OdorStart
    thisTrialSniffs(find(thisTrialSniffs(:,1)>=thisTrialOdorStart),6) = TrialInfo.Odor(trialID);
    % pull down any inhalations that happend after OdorEnd
    thisTrialSniffs(find(thisTrialSniffs(:,1)>=thisTrialStartStop(2)),6) = -1;

    SniffTimeStamps = vertcat(SniffTimeStamps, ...
        [thisTrialSniffs (trialID + 0*thisTrialSniffs(:,1))]);
end

end