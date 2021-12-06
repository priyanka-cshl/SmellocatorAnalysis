function [handles] = EventsPlotter(handles,BoxTag,RewardTag,TTLs,TuningTTLs)

for i = 1:3
    handles.([BoxTag,num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.([BoxTag,num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TTLs.(['Odor',num2str(i)])(:,1:2)';
    handles.([BoxTag,num2str(i),'Plot']).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat((10+str2double(handles.NumUnits.String))*[0 1 1 0]',size(ValveTS,2),1)];
    handles.([BoxTag,num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end

if ~isempty(TuningTTLs)
    % plot air periods for the tuning period
    i = 4;
    handles.([BoxTag,num2str(i),'Plot']) = fill(NaN,NaN,Plot_Colors(['Odor',num2str(i)]));
    hold on;
    handles.([BoxTag,num2str(i),'Plot']).EdgeColor = 'none';
    ValveTS = TuningTTLs(TuningTTLs(:,5)==1,[4 6])';
    handles.([BoxTag,num2str(i),'Plot']).Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat((10+str2double(handles.NumUnits.String))*[0 1 1 0]',size(ValveTS,2),1)];
    handles.([BoxTag,num2str(i),'Plot']).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
end

% plot Rewards
handles.(RewardTag) = plot(NaN, NaN, 'color',Plot_Colors('t'),'Linewidth',1.25);
tick_timestamps =  TTLs.Reward(:,1);
tick_x = [tick_timestamps'; tick_timestamps'; ...
    NaN(1,numel(tick_timestamps))]; % creates timestamp1 timestamp1 NaN timestamp2 timestamp2..
tick_x = tick_x(:);
tick_y = repmat( [0; 10+str2double(handles.NumUnits.String); NaN],...
    numel(tick_timestamps),1); % creates y1 y2 NaN y1 timestamp2..
set(handles.(RewardTag),'XData',tick_x,'YData',tick_y);