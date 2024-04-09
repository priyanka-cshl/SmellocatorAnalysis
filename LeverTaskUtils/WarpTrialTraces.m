function [Leverwarped,Odorwarped] = WarpTrialTraces(Traces, TrialInfo, whichTrials, tracebuffer)

if nargin<3
    whichTrials = 1:size(Traces.Lever,2); % all trials
end

% global SampleRate; % = 500; % samples/second
% global startoffset; % = 1; % seconds
% tracebuffer = SampleRate*0.5;

if nargin<4
    tracebuffer = 250; % 0.5 second
end

whichTraces = fieldnames(Traces);
warplength = 100; 

for t = 1:numel(whichTrials)
    thisTrial = whichTrials(t);
    t1 = find(diff(Traces.Trial{thisTrial})>0,1,'first');
    t2 = find(diff(Traces.Trial{thisTrial})<0,1,'last');
    tw = linspace(t1,t2,warplength);
    Leverwarped(t,:) = interp1(t1:t2, Traces.Lever{thisTrial}(t1:t2), tw);
    Odorwarped(t,:) = interp1(t1:t2, Traces.Motor{thisTrial}(t1:t2), tw);
end