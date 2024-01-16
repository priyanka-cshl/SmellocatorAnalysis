function [] = PlotSpikeAmplitudes_tester(varargin)

[Paths] = WhichComputer();

if ~isempty(varargin)
    if exist(varargin{1}) == 2
        handles.WhereSession.String = varargin{1};
    else
        MouseName = regexprep(varargin{1},'_(\w+)_processed.mat','');
        handles.WhereSession.String = fullfile(Paths.ProcessedSessions,MouseName,varargin{1});
    end
else
    handles.WhereSession.String = fullfile(Paths.ProcessedSessions,'O3/O3_20211005_r0_processed.mat');
end

MySession = handles.WhereSession.String;
% Load the relevant variables
load(MySession, 'TrialInfo', 'TTLs', 'ReplayTTLs', 'Tuning*', 'SingleUnits', 'TimestampAdjust', 'SniffTS', 'SniffTS_passive');

% concatenate SniffTS and SniffTS_passive and convert both to OEPS timebase
SniffTS = SniffTS + TimestampAdjust.ClosedLoop;
SniffTS = [SniffTS; (SniffTS_passive + TimestampAdjust.Passive)];

% load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
%     'startoffset', 'errorflags', 'SampleRate', ...
%     'TTLs', 'ReplayTTLs', 'Tuning*', 'SingleUnits');

handles.TrialInfo = TrialInfo;
handles.SingleUnits = SingleUnits;
handles.TimestampAdjuster = TimestampAdjust.ClosedLoop;

handles.spike_amplitudes.Value = 1;
whichUnit = 5;

% plot spike amplitudes
if handles.spike_amplitudes.Value
    subplot(2,1,1);
    thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
    %axes(handles.amplitudeaxes);
    hold off
    plot(handles.SingleUnits(whichUnit).spikes,thisunitamps,'.');
    hold on
    session_end = handles.TrialInfo.SessionTimestamps(end,2) + handles.TimestampAdjuster;
    line([session_end session_end],get(gca,'YLim'),'Color','k');
    
    subplot(2,1,2);
    % plot respiration phase
    thisunitspikes = handles.SingleUnits(whichUnit).spikes;
    for i = 1:numel(thisunitspikes)
        thisspike   = thisunitspikes(i);
        whichsniff  = find(SniffTS(:,1)<=thisspike,1,'last') ;
        if ~isempty(whichsniff)
            if thisspike <= SniffTS(whichsniff,3)
                handles.SingleUnits(whichUnit).sniffalignedspikes(i,1) = thisspike -  SniffTS(whichsniff,1); % latency
            else
                handles.SingleUnits(whichUnit).sniffalignedspikes(i,1) = NaN; 
            end
        else
            handles.SingleUnits(whichUnit).sniffalignedspikes(i,1) = NaN; 
        end
        
    end
    thisunitphases = handles.SingleUnits(whichUnit).sniffalignedspikes;
    hold off
    %plot(handles.SingleUnits(whichUnit).spikes,thisunitphases,'.');
    colormap(brewermap([10],'Dark2'));
    scatter(handles.SingleUnits(whichUnit).spikes,thisunitphases,12,thisunitamps/max(thisunitamps),'.');
    hold on
    line([session_end session_end],get(gca,'YLim'),'Color','k');
end