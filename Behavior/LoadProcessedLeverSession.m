
function [TracesOut, whichTraces, SegmentedLever, thisSniffParams, TrialTimeStamps, TrialIndices, startspositive, TrialInfo, SampleRate, TargetZones] = ...
            LoadProcessedLeverSession(WhereSession)

% WhereSession = '/home/priyanka/Dropbox/Smellocator_Behavior_2/Q5/Q5_20221028_r0_processed.mat';
% WhereSession = '/home/priyanka/Dropbox/Smellocator_Behavior_2/Q3/Q3_20221012_r0_processed.mat';
% WhereSession = '/home/priyanka/Dropbox/Smellocator_Behavior_2/Q4/Q4_20221129_r0_processed.mat';
% WhereSession = '/home/priyanka/Dropbox/Smellocator_Behavior_2/Q9/Q9_20221108_r0_processed.mat';


%%
load(WhereSession, 'Traces', 'TrialInfo', 'TargetZones','startoffset', 'SampleRate');
[TracesOut, whichTraces] = ConcatenateTraces2Matrix(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);

% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(whichTraces,'Timestamps')));

%% Calculate Trial On-Off timestamps
TrialColumn = TracesOut(:,find(strcmp(whichTraces,'TrialState')));
TrialColumn(TrialColumn~=0) = 1; % make logical
TrialOn = find(diff([0; TrialColumn])>0);
TrialOff =  find(diff(TrialColumn)<0)+1;

% account for cases where acquisition started/ended in between a trial
while TrialOn(1)>TrialOff(1)
    TrialOff(1,:) = [];
end
while TrialOn(end)>TrialOff(end)
    TrialOn(end,:) = [];
end

TrialIndices = [TrialOn TrialOff];
TrialTimeStamps = [Timestamps(TrialOn,1) Timestamps(TrialOff,1)];

startspositive = median(TracesOut(TrialIndices(:,1)+1,5))>0;

if size(TrialTimeStamps,1) == size(TrialInfo.OdorStart,1)
    
    TrialTimeStamps(:,1) = TrialTimeStamps(:,1) + TrialInfo.OdorStart(:,1);
    TrialIndices(:,3)    = TrialIndices(:,1) + round(SampleRate*TrialInfo.OdorStart(:,1));
end

%% Sniffing specifc
% add a filtered sniff trace
TracesOut(:,9)     = FilterThermistor(TracesOut(:,2));
whichTraces{9,1}   = 'SniffsFiltered';

% make a digital sniff trace
load(WhereSession,'CuratedSniffTimestamps');
if exist('CuratedSniffTimestamps','var')

    if size(CuratedSniffTimestamps,2) < 10
        CuratedSniffTimestamps(:,10) = 0;
    end
    LocationSniffs  = TracesOut(:,9)*nan;
    DigitalSniffs   = TracesOut(:,9)*0;
    SegmentedLever  = zeros(size(CuratedSniffTimestamps,1)-1,8);
    for n = 1:size(CuratedSniffTimestamps,1)-1
        idx = CuratedSniffTimestamps(n,8:9); % trace indices of inhalation start and end
        if CuratedSniffTimestamps(n,10) == -1
            DigitalSniffs(idx(1):idx(2)) = -1;
        else
            DigitalSniffs(idx(1):idx(2)) = 1;
        end
        if ~any(isnan(CuratedSniffTimestamps(:,4)))
            location = CuratedSniffTimestamps(n,4);
            LocationSniffs(idx(1):idx(2)) = location;
        end
        SegmentedLever(n,1:4) = [CuratedSniffTimestamps(n,1:2) TracesOut(idx,1)']; % ts and lever values at inhalation start and end
        SegmentedLever(n,5:6) = [mean(CuratedSniffTimestamps(n,1:2)) mean(TracesOut(idx(1):idx(2),1))]; % mean ts and lever value at inhalation
        idx(3) = CuratedSniffTimestamps(n+1,8); % next sniff start
        SegmentedLever(n,7) = mean([CuratedSniffTimestamps(n,2) CuratedSniffTimestamps(n+1,1)]) ;
        SegmentedLever(n,8) = mean(TracesOut(idx(2):idx(3),1)); % mean ts and lever value at 'rest of sniff'
    end

    TracesOut(:,10)     = DigitalSniffs;
    whichTraces{10,1}   = 'SniffsDigitized';
    TracesOut(:,11)     = LocationSniffs;
    whichTraces{11,1}   = 'SniffsLocationed';
    TracesOut(:,12)     = TracesOut(:,1); % Lever
    TracesOut(find(~TracesOut(:,10)),12) = nan;
    whichTraces{12,1}   = 'LeverGated';
else
    SegmentedLever = [];
end

%% plotting sniff segmented lever moves
whichsniffs = find(CuratedSniffTimestamps(:,6)>0);
whichsniffs(end,:) = [];

for s = 1:numel(whichsniffs)
    thisSniffLocation   = CuratedSniffTimestamps(whichsniffs(s),4);
    prevSniffLocation   = CuratedSniffTimestamps(whichsniffs(s)-1,4);
    thisSniffLever      = [ SegmentedLever(whichsniffs(s),4)-SegmentedLever(whichsniffs(s),3), ...
                            SegmentedLever(whichsniffs(s)+1,3)-SegmentedLever(whichsniffs(s),4), ...
                            SegmentedLever(whichsniffs(s)+1,3)-SegmentedLever(whichsniffs(s),3) ];
    thisSniffTimes      = [ CuratedSniffTimestamps(whichsniffs(s),2) - CuratedSniffTimestamps(whichsniffs(s),1), ...
                            CuratedSniffTimestamps(whichsniffs(s),3) - CuratedSniffTimestamps(whichsniffs(s),2), ...
                            CuratedSniffTimestamps(whichsniffs(s),3) - CuratedSniffTimestamps(whichsniffs(s),1) ];
    thisSniffAngles     = atand(thisSniffLever./thisSniffTimes);
    thisSniffTrialMid = mean(TrialTimeStamps(CuratedSniffTimestamps(whichsniffs(s),7),1:2));
    thisSniffParams(s,:)= [thisSniffLocation thisSniffLever thisSniffTimes thisSniffAngles thisSniffTrialMid prevSniffLocation];
    % odorloc 
    % leverdisinh leverdisexh leverdissnf durationinh durationexh durationsnf
    % leveranginh leverangexh leverangsnf trialstart  odorlocprevious
end



end
%% for every trial
% figure, for i = 1:3; subplot(1,3,i); scatter(thisSniffParams(:,1),thisSniffParams(:,i+1+6)); end
% for t = 1:size(TrialIndices,1)
% 
% end


%% plot odor boxes on the behavior plot
% axes(handles.BehaviorPlot);
% hold off
% for i = 1:4
%     handles.(['Trial',num2str(i),'Plot']) = fill(NaN,NaN,handles.plotcolors.(['Odor',num2str(i)]));
%     hold on;
%     handles.(['Trial',num2str(i),'Plot']).EdgeColor = 'none';
%     ValveTS = TrialTimeStamps((TrialInfo.Odor==i),1:2)';
%     if ~isempty(ValveTS)
%         handles.(['Trial',num2str(i),'Plot']).Vertices = [ ...
%             reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
%             repmat([0 5 5 0]',size(ValveTS,2),1)];
%         handles.(['Trial',num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
%     end
% end