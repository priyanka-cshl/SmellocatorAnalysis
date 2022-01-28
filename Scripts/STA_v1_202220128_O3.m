MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up
%MySession = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/O3/O3_20211005_r0_processed.mat';
%MySession = '/mnt/data/Processed/Behavior/PCX4/PCX4_20210721_r0_processed.mat'; % session path - leave empty to get browser pop up
LoadProcessedSession; % loads relevant variables

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

N = size(SingleUnits,2);
% sort units by tetrode - to match session viewer
clear foo
for i = 1:N
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id];
end
[~,SortedByTetrodes] = sort(foo(:,1));
UnitOrder = SortedByTetrodes;

% SampleRate = behavior sample rate;
[TracesOut] = ConcatenateTraces(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);
% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

ChosenUnits = [58 35 34 55 21];
STAwindow = [-1 1];
snippetlength = diff(STAwindow)*SampleRate;
F = fieldnames(TracesOut);
% for every unit
for i = 1:numel(ChosenUnits)
    MyUnit = ChosenUnits(i);
    % for every spike (in behavior timebase
        thisunitspiketimes = SingleUnits(MyUnit).spikes - TimestampAdjuster + STAwindow(1);
        % ignore spikes that precede behavior start
        thisunitspiketimes(thisunitspiketimes<0) = [];
        % ignore spikes that follow ebhavior stop
        thisunitspiketimes(thisunitspiketimes>=(Timestamps(end)-STAwindow(2))) = [];
        window = [];
        Snippets = [];
        for j = 1:numel(thisunitspiketimes)
            idx1 = find(Timestamps>=thisunitspiketimes(j),1,'first');
            idx2 = idx1 + snippetlength;
            for k = 1:numel(F)-1
%                Snippets(j,:,k) = TracesOut.(F{k}){1}(idx1:idx2);
                mysnippet = TracesOut.(F{k}){1}(idx1:idx2);
                if j == 1
                    Snippets(:,k) = mysnippet;
                else
                    Snippets(:,k) = Snippets(:,k) + mysnippet;
                end
            end
        end
        % get the average trace
        STA(:,:,i) = Snippets/j;
%         for k = 1:numel(F)-1
%             STA(:,k,i) = Snippets
%         end
end
    
% TraceNames = {'Lever' 'Motor' 'Sniffs' 'Licks' 'Rewards' 'Trial' 'TargetZone'};

