function [TracesOut] = ConcatenateTraces(Traces, whichTrials, traceOverlap)

% global SampleRate; % = 500; % samples/second
% global startoffset; % = 1; % seconds
% traceOverlap = SampleRate*startoffset;

whichTraces = fieldnames(Traces);

% get all traces and concatenate them
    for j = 1:size(whichTraces,1)
        temp = cellfun(@(x) ...
            x(1:end-traceOverlap), Traces.(whichTraces{j})(whichTrials), ...
            'UniformOutput', false);
        
        % add in the overlap for the very last trial if needed
        % also make sure that there are atleast startoffset*SampleRate
        % samples after TrialOFF
        TracesOut.(whichTraces{j}) = {[cell2mat(temp(:)); ...
            Traces.(whichTraces{j}){whichTrials(end)}(end-traceOverlap+1:end,1)]};
    end
end