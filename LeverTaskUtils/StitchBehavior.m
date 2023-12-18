function [TracesOut, TrialInfoOut, SniffTSOut] = StitchBehavior(TracesIn,TrialInfoIn,SniffTSIn,OEPSIndexes)

TracesOut       = TracesIn{1};
TrialInfoOut    = TrialInfoIn{1};
SniffTSOut      = SniffTSIn{1};
trialselapsed   = size(TrialInfoIn{1}.TrialID,2); % trials in session 1
whichTraces     = fieldnames(TracesIn{1});
trialfields     = fieldnames(TrialInfoIn{1});

TrialInfoOut.OepsID     = (OEPSIndexes(1,1):OEPSIndexes(1,2))';
for n = 2:size(TracesIn,2) % no. of sessions
    % how much time elapsed between the two sessions - no. to be added to
    % adjust timestamps for the new session w.r.t to Session I 
    session_gap = TrialInfoIn{1}.SessionTimestamps(1,2) + OEPSIndexes(n,3) - OEPSIndexes(1,3) - TrialInfoIn{n}.SessionTimestamps(1,2);
    
    % Traces - concatenate and adjust SessionTimetamps
    for j = 1:size(whichTraces,1)
        if strcmp(whichTraces{j},'Timestamps')
            TracesOut.(whichTraces{j}) = [TracesOut.(whichTraces{j}) cellfun(@(x) session_gap+x, TracesIn{n}.Timestamps, 'UniformOutput', false)];
        else
            TracesOut.(whichTraces{j}) = [TracesOut.(whichTraces{j}) TracesIn{n}.(whichTraces{j})];
        end
    end
    
    % TrialInfo - concatenate and adjust SessionTimetamps
    for k = 1:size(trialfields,1)
        if ~iscell(TrialInfoIn{1}.(trialfields{k}))
            % find whether to append along rows or columns
            if size(TrialInfoIn{1}.(trialfields{k}),1) == trialselapsed
                if strcmp(trialfields{k},'SessionTimestamps')
                    TrialInfoOut.(trialfields{k}) = vertcat(TrialInfoOut.(trialfields{k}), session_gap + TrialInfoIn{n}.(trialfields{k}));
                else
                    TrialInfoOut.(trialfields{k}) = vertcat(TrialInfoOut.(trialfields{k}), TrialInfoIn{n}.(trialfields{k}));
                end
            else
                TrialInfoOut.(trialfields{k}) = horzcat(TrialInfoOut.(trialfields{k}), TrialInfoIn{n}.(trialfields{k}));
            end
        else
            if size(TrialInfoIn{1}.(trialfields{k}),1) == trialselapsed
                TrialInfoOut.(trialfields{k}) = [TrialInfoOut.(trialfields{k}); TrialInfoIn{n}.(trialfields{k})];
            else
                TrialInfoOut.(trialfields{k}) = [TrialInfoOut.(trialfields{k}) TrialInfoIn{n}.(trialfields{k})];
            end
        end
        
    end
    % create a vector for corresponding OEPS trial IDs
    TrialInfoOut.OepsID = vertcat(TrialInfoOut.OepsID, (OEPSIndexes(n,1):OEPSIndexes(n,2))');
    
    % SniffTimestamps - concatenate and adjust SessionTimetamps
    SniffTSOut = vertcat(SniffTSOut, [(session_gap + SniffTSIn{n}(:,1:3)) SniffTSIn{n}(:,4)]);
    
end

TrialInfoOut.SessionSplits = OEPSIndexes(:,1:2);

end