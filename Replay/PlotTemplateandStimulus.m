
function [MyTraces,timestamps,PSTH,Raster] = PlotTemplateandStimulus(Replay, TrialInfo)

global SampleRate;
global MyFileName;
global TargetZones;
global startoffset;


allreplays = 1:numel(Replay.TemplateTraces.TrialIDs);
PassiveReplays = Replay.TTLs.TrialID(Replay.TTLs.TrialID>TrialInfo.TrialID(end));

whichreplay = 1;
    
%% 1. Behavior

% get trace length from the replay traces
tracelength = size(Replay.ReplayTraces.Lever{whichreplay},1);

TraceNames = {'Lever' 'Motor' 'Sniffs' 'Licks' 'Rewards' 'Trial' 'TargetZone'};
for i = 1:numel(TraceNames)
    temptrace = Replay.TemplateTraces.(TraceNames{i}){whichreplay};
    temptrace = temptrace(~isnan(temptrace));
    MyTraces(:,i,1) = temptrace(1:tracelength,1);
end
    
% common traces for both template and replays
% Trial, TargetZone and Timestamps

Trial = MyTraces(:,6,1);
% Trial column has -ve values that indicate odorON periods
% ignore them for plotting
%Trial(Trial<0) = 0;
Trial = abs(Trial);

TZ_Up   = MyTraces(:,7,1);
TZ_Low  = TZ_Up;
% change TZ vector to TZ_lims
MyZones = flipud(unique(TZ_Up(TZ_Up~=0)));
for i = 1:numel(MyZones)
    whichzone = MyZones(i);
    TZ_Up(TZ_Up==whichzone) = TargetZones(find(TargetZones(:,2)==whichzone),1);
    TZ_Low(TZ_Low==whichzone) = TargetZones(find(TargetZones(:,2)==whichzone),3);
end
TZ = [TZ_Low TZ_Up];
MotorTZ = 0*TZ;
MotorTZ(:,1) = -8+MotorTZ(:,1);
MotorTZ(:,2) = 8+MotorTZ(:,2);

timestamps = (1:tracelength)'/SampleRate;
    
% Plotting the behavior
H1 = figure; % Lever traces
figure(H1);

subplot(2,1,1);
PlotBehavior(timestamps,MyTraces(:,1,1),MyTraces(:,3,1),MyTraces(:,4,1),MyTraces(:,5,1),Trial,TZ);
set(gca,'YLim',[-0.4 8],'YTick',[0 5],'TickDir','out','XLim',[0 round(timestamps(end))],...
    'XTick',[]);

subplot(2,1,2);
imagesc(timestamps,1,(-MyTraces(:,2,1))');
colormap(brewermap([100],'*rdbu')); set(gca,'CLim',160*[-1 1]);
set(gca,'YLim',[0 2],'YTick',[],'TickDir','out','XLim',[0 round(timestamps(end))],...
    'XTick',[]);

% overlay Trial OFF periods
hold on
handles.trial_off = fill(NaN,NaN,[1 1 1]);
handles.trial_off.EdgeColor = 'none';
[handles] = PlotToPatch_TrialOFFhack(handles, Trial, timestamps, [-1 2],1);

end