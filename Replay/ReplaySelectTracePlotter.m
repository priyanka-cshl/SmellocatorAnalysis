function [] = ReplaySelectTracePlotter(TemplateTraces, ActiveReplayTraces, PassiveReplayTraces, TraceTag)

if nargin < 4
    TraceTag = 'Sniffs';
end

n(1) = size(TemplateTraces.Timestamps{1},2); % Template replays
n(2) = size(ActiveReplayTraces.Timestamps{1},2); % Active replays
n(3) = size(PassiveReplayTraces.Timestamps,2); % Passive replays

figure;
% plot template
for i = 1:n(1)
    subplot(sum(n),1,i)
    plot(TemplateTraces.Timestamps{1}(10:end,i) - ...
        TemplateTraces.Timestamps{1}(10,i) + .02, ...
        TemplateTraces.(TraceTag){1}(10:end,i), 'color', Plot_Colors('k'));
    set(gca,'XAxisLocation','top');
end
% plot active replays
for i = 1:n(2)
    subplot(sum(n),1,i+n(1))
    plot(ActiveReplayTraces.Timestamps{1}(10:end,i) - ...
        ActiveReplayTraces.Timestamps{1}(10,i) + .02, ...
        ActiveReplayTraces.(TraceTag){1}(10:end,i), 'color', Plot_Colors('t'));
    set(gca,'XTick',[]);
end
% plot passive replays
for i = 1:n(3)
    subplot(sum(n),1,i+n(1)+n(2))
    plot(PassiveReplayTraces.Timestamps{i} - ...
        PassiveReplayTraces.Timestamps{i}(1), ...
        PassiveReplayTraces.(TraceTag){i}, 'color', Plot_Colors('r'));
    set(gca,'XTick',[]);
end
end