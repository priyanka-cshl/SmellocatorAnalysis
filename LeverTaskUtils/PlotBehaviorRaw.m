function [] = PlotBehaviorRaw(MyFilePath)

[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
timestamps  = MyData(:,1);
Lever       = MyData(:,4);
Sniffs      = MyData(:,15);
Rewards     = MyData(:,9);
Licks       = MyData(:,10);
Trial       = MyData(:,6);
TZ          = MyData(:,[2 3]);
OdorBoxHeight = 5;

global SampleRate;
SampleRate = 500;
global startoffset;
startoffset = 1; % in seconds
global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips]
errorflags = [0 0 0 0];

[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags);
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

%% Plots
% Trials
handles.trial_on_1 = fill(NaN,NaN,[.8 .8 .8]);
hold on;
handles.trial_on_1.EdgeColor = 'none';
handles.trial_on_2 = fill(NaN,NaN,[0.8941    0.9412    0.9020]);
%handles.trial_on_2 = fill(NaN,NaN,[.8 .8 .8]);
handles.trial_on_2.EdgeColor = 'none';
handles.trial_on_3 = fill(NaN,NaN,[0.8706    0.9216    0.9804]);
%handles.trial_on_3 = fill(NaN,NaN,[.8 .8 .8]);
handles.trial_on_3.EdgeColor = 'none';
handles.trial_on_4 = fill(NaN,NaN,[0.93    0.84    0.84]);
handles.trial_on_4.EdgeColor = 'none';

% for plotting Trial as TF
handles.trial_tf = fill(NaN,NaN,[1 1 1]);
handles.trial_tf.EdgeColor = 'none';
colormap(brewermap([100],'rdbu'));
set(gca,'CLim',120*[-1 1]);

% whiten trial off periods
handles.trial_off = fill(NaN,NaN,[1 1 1]);
handles.trial_off.EdgeColor = 'none';

if ~isempty(Trial)
    if ~isempty(TZ)
        [handles] = PlotToPatch_TrialTF(handles, Trial, timestamps, [0 OdorBoxHeight],mean(TZ,2));
        [handles] = PlotToPatch_Trial(handles, Trial, timestamps, [0 -1],1);
    else
        [handles] = PlotToPatch_Trial(handles, Trial, timestamps, [0 OdorBoxHeight],1);
    end
    
    % old version
    % [handles] = PlotToPatch_Trial(handles, Trial, timestamps, [0 OdorBoxHeight],1);
    % [handles] = PlotToPatch_TrialOFFhack(handles, Trial, timestamps, [0 OdorBoxHeight],1);
end

% TargetZone
if ~isempty(TZ)
    %handles.targetzone = fill(NaN,NaN,[1 1 0],'FaceAlpha',0.2);
    handles.targetzone = fill(NaN,NaN,[1 1 1]);
    handles.targetzone.EdgeColor = 'none';

    [handles.targetzone] = PlotToPatch_TargetZone(handles.targetzone, abs(TZ), timestamps);
end

% Lever
if ~isempty(Lever)
    handles.lever_plot = plot(NaN, NaN,'k','Linewidth',1);
    set(handles.lever_plot,'XData',timestamps,'YData',Lever);
end

% Sniffs
if ~isempty(Sniffs)
    handles.respiration_plot = plot(NaN, NaN,'color',Plot_Colors('r'));
    set(handles.respiration_plot,'XData',timestamps,'YData',6.5+ 2*(Sniffs/max(Sniffs)));
end

% Rewards
if ~isempty(Rewards)
    handles.reward_plot = plot(NaN, NaN, 'color',Plot_Colors('t'),'Linewidth',1.25);
    tick_timestamps =  timestamps(Rewards==1);
    tick_x = [tick_timestamps'; tick_timestamps'; ...
        NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
    tick_x = tick_x(:);
    tick_y = repmat( [0; 5; NaN],...
        numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
    set(handles.reward_plot,'XData',tick_x,'YData',tick_y);
end

% Licks
if ~isempty(Licks)
    handles.lick_plot = plot(NaN, NaN, 'color',Plot_Colors('o'),'Linewidth',1); %licks
    tick_timestamps =  timestamps(Licks==1);
    tick_x = [tick_timestamps'; tick_timestamps'; ...
        NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
    tick_x = tick_x(:);
    tick_y = repmat( [5.5; 6; NaN],...
        numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
    set(handles.lick_plot,'XData',tick_x,'YData',tick_y);
end

% plot OdorStart
plot(TrialInfo.OdorStart(:,2)+TrialInfo.SessionTimestamps(:,1),-0.25,'.k');
end
