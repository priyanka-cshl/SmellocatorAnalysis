function [timestamps,Lever,Sniffs,Licks,Rewards,Trial,TZ,OdorBoxHeight] = ...
    PrepBehaviorForPlotting(MyData, DataTags)

timestamps  = MyData(:,find(cellfun(@isempty,regexp(DataTags,'Timestamps'))==0));
Lever       = MyData(:,find(cellfun(@isempty,regexp(DataTags,'Lever'))==0));
Sniffs      = MyData(:,find(cellfun(@isempty,regexp(DataTags,'respiration'))==0));
Licks       = MyData(:,find(cellfun(@isempty,regexp(DataTags,'Licks'))==0));
Rewards     = MyData(:,find(cellfun(@isempty,regexp(DataTags,'Rewards'))==0));
Trial       = MyData(:,find(cellfun(@isempty,regexp(DataTags,'TrialON'))==0));
TZ(:,1)     = MyData(:,find(cellfun(@isempty,regexp(DataTags,'TZoneHighLim'))==0));
TZ(:,2)     = MyData(:,find(cellfun(@isempty,regexp(DataTags,'TZoneLowLim'))==0));
MotorTZ = 0*TZ;
MotorTZ(:,1) = -8+MotorTZ(:,1);
MotorTZ(:,2) = 8+MotorTZ(:,2);

end
