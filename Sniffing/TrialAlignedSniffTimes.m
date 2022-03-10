function [AlignedSniffs] = TrialAlignedSniffTimes(SniffTimes,Trial,SampleRate)
if size(Trial,2) == 1
    % detect Trial On and Off
    Trial(Trial>0) = 1;
    Trial = vertcat(0,Trial,0);
    TS = [1+find(diff(Trial)==1) find(diff(Trial)==-1)];
else
    TS = Trial; % user already input trial On-Off times in indices
end

for i = 1:size(TS,1)
    x1 = TS(i,1) - SampleRate;
    x2 = TS(i,2) + SampleRate;
    thisTrialSniffs = intersect(find(SniffTimes(:,1)>x1),...
                                find(SniffTimes(:,1)<=x2));
                            
    AlignedSniffs{i,1} = {(SniffTimes(thisTrialSniffs,1) - TS(i,1))/SampleRate};
    if size(SniffTimes,2)==2
        AlignedSniffs{i,2} = {(SniffTimes(thisTrialSniffs,2) - TS(i,1))/SampleRate};
    end
end

end