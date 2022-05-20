% recreate session plot in session-plotting style
% three different trial plots - for different odor
function [MyFig,MyHandle] = RecreateSession_v2(MyData)

if ~isempty(regexp(MyData,'.mat','match'))
    [MyData] = ReadSessionData(MyData);
end

MyFig = figure('Name','Session','NumberTitle','off');

%% initialize plots
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

% handles.fake_target_plot = plot(NaN, NaN, 'color',[.7 .7 .7]);
handles.targetzone = fill(NaN,NaN,[1 1 0],'FaceAlpha',0.2);
handles.targetzone.EdgeColor = 'none';
handles.faketargetzone = fill(NaN,NaN,[1 0 1],'FaceAlpha',0.2);
handles.faketargetzone.EdgeColor = 'none';

% for plotting Trial as TF
handles.trial_tf = fill(NaN,NaN,[1 1 1]);
handles.trial_tf.EdgeColor = 'none';
colormap(brewermap([100],'rdbu'));
set(gca,'CLim',120*[-1 1]);

% whiten trial off periods
handles.trial_off = fill(NaN,NaN,[1 1 1]);
handles.trial_off.EdgeColor = 'none';

handles.lever_DAC_plot = plot(NaN, NaN,'k','Linewidth',0.75); %lever rescaled
handles.Camera_plot = plot(NaN, NaN,'k');
handles.respiration_plot = plot(NaN, NaN,'color',Plot_Colors('r')); %lever rescaled
handles.stimulus_plot = plot(NaN, NaN, 'color',Plot_Colors('r'),'Linewidth',1); % target odor location
handles.reward_plot = plot(NaN, NaN, 'color',Plot_Colors('t'),'Linewidth',1.25); %rewards
handles.lick_plot = plot(NaN, NaN, 'color',Plot_Colors('o'),'Linewidth',1); %licks

set(gca,'YLim',[-0.4 8]);

MyHandle = get(gca);
%% Update plots

% lever positions, motor locations 
% [a,b] = butter(3,100/250,'low');
% mylever = filter(a,b,MyData(:,4));
% set(handles.lever_DAC_plot,'XData',MyData(:,1),'YData',mylever);
set(handles.lever_DAC_plot,'XData',MyData(:,1),'YData',MyData(:,4));

set(handles.stimulus_plot,'XData',MyData(:,1),'YData',...
   1*(MyData(:,5)- 0) );
if size(MyData,2)>=15
    set(handles.respiration_plot,'XData',MyData(:,1),'YData',6.5+ 2*(MyData(:,15)/max(MyData(:,15))));
end
if size(MyData,2)>15
    set(handles.Camera_plot,'XData',MyData(:,1),'YData',7+ 0.5*MyData(:,16));
end
% trial_on
[handles] = PlotToPatch_Trial(handles, MyData(:,6), MyData(:,1), [0 5],1);
[handles] = PlotToPatch_TrialOFFhack(handles, MyData(:,6), MyData(:,1), [0 5],1);
[handles.targetzone] = PlotToPatch_TargetZone(handles.targetzone, MyData(:,2:3), MyData(:,1));
[handles.faketargetzone] = PlotToPatch_TargetZone(handles.faketargetzone, MyData(:,11:12), MyData(:,1));

% % in_target_zone, in_reward_zone
% [handles.in_target_zone_plot] = PlotToPatch(handles.in_target_zone_plot, MyData(:,7), MyData(:,1), [-1 0],1);
% [handles.in_reward_zone_plot] = PlotToPatch(handles.in_reward_zone_plot,  MyData(:,8), MyData(:,1), [-1 -0.2],1);

% rewards
tick_timestamps =  MyData(MyData(:,9)==1,1);
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; 5; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.reward_plot,'XData',tick_x,'YData',tick_y);

% licks
tick_timestamps = MyData(MyData(:,10)==1,1);
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [5.5; 6; NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.lick_plot,'XData',tick_x,'YData',tick_y);

% set(handles.fake_target_plot,'XData',TotalTime(indices_to_plot),'YData',...
%         handles.PerturbationSettings.Data(4) + 0*TotalTime(indices_to_plot));

