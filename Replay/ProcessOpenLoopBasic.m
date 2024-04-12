function [] = ProcessOpenLoopBasic(Replay, SampleRate, TargetZones, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('whichreplays', [], @(x) isnumeric(x));
params.addParameter('plotbehavior', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plottrials', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('ShadeErrorbars', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('TrialHeight', [0 50], @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
allreplays = params.Results.whichreplays;
plotreplayfigs = params.Results.plotbehavior;
plottrialboxes = params.Results.plottrials;
Ylims = params.Results.TrialHeight;

%% just plotting replay behavior or the just the trial boxes - no ephys
whichreplay = 1;

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
Trial(Trial<0) = 0;

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

trials_per_replay = size(Replay.ReplayTraces.TrialIDs{whichreplay},1);

% data (specific to open loop)
for i = 1:numel(TraceNames)-2
    MyTraces(:,i,1+(1:trials_per_replay)) = Replay.ReplayTraces.(TraceNames{i}){whichreplay};
end

% Plotting the behavior
if plotreplayfigs
    H1 = figure; % Lever traces
    H2 = figure; % Motor - odor location - traces
    
    for j = 1:(trials_per_replay+1)
        
        figure(H1);
        subplot(trials_per_replay+1,1,j);
        PlotBehavior(timestamps,MyTraces(:,1,j),MyTraces(:,3,j),MyTraces(:,4,j),MyTraces(:,5,j),Trial,TZ);
        if j > trials_per_replay
            set(gca,'YLim',[-0.4 8],'YTick',[0 5],'TickDir','out','XLim',[0 round(timestamps(end))]);
        else
            set(gca,'YLim',[-0.4 8],'YTick',[0 5],'TickDir','out','XLim',[0 round(timestamps(end))],...
                'XTick',[]);
        end
        figure(H2);
        subplot(trials_per_replay+1,1,j);
        PlotBehavior(timestamps,(MyTraces(:,2,j)+100)/40,[],[],[],Trial,(MotorTZ+100)/40);
        if j > trials_per_replay
            set(gca,'YLim',[-0.4 5.4],'YTick',[0 5],'TickDir','out','XLim',[0 round(timestamps(end))]);
        else
            set(gca,'YLim',[-0.4 5.4],'YTick',[0 5],'TickDir','out','XLim',[0 round(timestamps(end))],...
                'XTick',[]);
        end
    end
end

% plotting just trial boxes for physiology
if plottrialboxes
    % plot the trial structure
    PlotBehavior(timestamps,[],[],[],[],Trial,[],Ylims(2));
    set(gca,'YLim',Ylims,'YTick',Ylims,'TickDir','out','XLim',[0 round(timestamps(end))], ...
        'XTick', []);
    axis manual;
end
