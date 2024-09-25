
function [MyTraces,timestamps,PSTH,Raster] = ProcessOpenLoopTrials_v2(Replay, TrialInfo, SingleUnits, TTLs, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('whichreplays', [], @(x) isnumeric(x));
params.addParameter('plotfigures', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotephys', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('savefigures', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('whichunits', [], @(x) isnumeric(x));
params.addParameter('PSTHsmooth', 100, @(x) isnumeric(x));
params.addParameter('UnitsPerFig', 8, @(x) isnumeric(x));
params.addParameter('PlotOpenLoop', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('PlotPassive', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('ShadeErrorbars', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
allreplays = params.Results.whichreplays;
plotreplayfigs = params.Results.plotfigures;
savereplayfigs = params.Results.savefigures;
plotephysfigs = params.Results.plotephys;
whichUnits = params.Results.whichunits;
smoothingfactor = params.Results.PSTHsmooth;
units_per_fig = params.Results.UnitsPerFig;
plotOL = params.Results.PlotOpenLoop;
plotPR = params.Results.PlotPassive;
ShadeErrorbars = params.Results.ShadeErrorbars;

global SampleRate;
global MyFileName;
global TargetZones;
global startoffset;

if isempty(allreplays)
    allreplays = 1:numel(Replay.TemplateTraces.TrialIDs);
end

if ~isempty(SingleUnits) && isempty(whichUnits)
    whichUnits = 1:size(SingleUnits,2);
end

TraceTag = 'PassiveReplayTraces';
PassiveReplays = Replay.TTLs.TrialID(Replay.TTLs.TrialID>TrialInfo.TrialID(end));

%     %% 1. Behavior
%
%     % get trace length from the replay traces
%     tracelength = size(Replay.(TraceTag).Lever{whichreplay},1);
%
%     TraceNames = {'Lever' 'Motor' 'Sniffs' 'Licks' 'Rewards' 'Trial' 'TargetZone'};
%     for i = 1:numel(TraceNames)
%         temptrace = Replay.TemplateTraces.(TraceNames{i}){whichreplay};
%         temptrace = temptrace(~isnan(temptrace));
%         MyTraces(:,i,1) = temptrace(1:tracelength,1);
%     end
%
%     % common traces for both template and replays
%     % Trial, TargetZone and Timestamps
%
%     Trial = MyTraces(:,6,1);
%     % Trial column has -ve values that indicate odorON periods
%     % ignore them for plotting
%     Trial(Trial<0) = 0;
%
%     TZ_Up   = MyTraces(:,7,1);
%     TZ_Low  = TZ_Up;
%     % change TZ vector to TZ_lims
%     MyZones = flipud(unique(TZ_Up(TZ_Up~=0)));
%     for i = 1:numel(MyZones)
%         whichzone = MyZones(i);
%         TZ_Up(TZ_Up==whichzone) = TargetZones(find(TargetZones(:,2)==whichzone),1);
%         TZ_Low(TZ_Low==whichzone) = TargetZones(find(TargetZones(:,2)==whichzone),3);
%     end
%     TZ = [TZ_Low TZ_Up];
%     MotorTZ = 0*TZ;
%     MotorTZ(:,1) = -8+MotorTZ(:,1);
%     MotorTZ(:,2) = 8+MotorTZ(:,2);
%
%     timestamps = (1:tracelength)'/SampleRate;
%
%     trials_per_replay = size(Replay.(TraceTag).TrialIDs{whichreplay},1);
%
%     % data (specific to open loop)
%     for i = 1:numel(TraceNames)-2
%         MyTraces(:,i,1+(1:trials_per_replay)) = Replay.(TraceTag).(TraceNames{i}){whichreplay};
%     end

%% Ephys
PSTH = []; % dimensions: [trial x samples x units]
if ~isempty(SingleUnits)
    %units_per_fig = 5;
    if plotreplayfigs || savereplayfigs
        figure;
    end

    for i = 1:numel(whichUnits) % for every cell
        MyUnit = whichUnits(i);
        allspikes = SingleUnits(MyUnit).spikes;

        FRmax = 25; % for plotting

        for x = 1:numel(allreplays) % for every unique replay stretch
            whichreplay = allreplays(x);

            % get trial IDs for the closed loop template trials
            TemplateTrials = Replay.TemplateTraces.TrialIDs{whichreplay};
            nsubtrials = numel(TemplateTrials);
            ReplayTrialLength = size(Replay.(TraceTag).Lever{whichreplay},1)/SampleRate; % in seconds
            TemplateTrialLength = size(Replay.TemplateTraces.Lever{whichreplay},1)/SampleRate; % in seconds
            FirstTrialDuration = TrialInfo.Duration(TemplateTrials(1));

            if plotephysfigs || savereplayfigs
                if mod(i,units_per_fig)
                    FRplot = 2*mod(i,units_per_fig);
                else
                    FRplot = 2*units_per_fig;
                end
                Rasterplot = FRplot - 1;

                subplot(units_per_fig,2,Rasterplot);
                % plot the trial structure
                %PlotBehavior(timestamps,[],[],[],[],Trial,[],numel(PassiveReplays)+trials_per_replay+1);
                tracelength = size(Replay.(TraceTag).Lever{whichreplay},1);
                timestamps = (1:tracelength)'/SampleRate;
                set(gca,'YLim',...
                    [-0.4 (2*numel(allreplays))+0.4],...
                    'YTick',[0 5],'TickDir','out','XLim',[0 round(timestamps(end))], ...
                    'XTick', []);
                axis manual;
                hold on;

%                 subplot(units_per_fig,2,FRplot);
%                 % plot the trial structure
%                 %PlotBehavior(timestamps,[],[],[],[],Trial,[],FRmax);
%                 set(gca,'YLim',[0 FRmax],'YTick',[0 50 100],'TickDir','out','XLim',[0 round(timestamps(end))], ...
%                     'XTick', []);
%                 axis manual;
%                 hold on;
            end

            % Get Spikes
            thisTrialSpikeTimes = [];
            % first collate spikes for the original close loop stretch of
            % template trials
            thisTrialSpikeTimes = allspikes(ismember(SingleUnits(MyUnit).trialtags,TemplateTrials));
            % append in the beginning the spikes in the startoffset preceeding trial 1
            tstart = TTLs.Trial(TemplateTrials(1),1);
            thisTrialSpikeTimes = vertcat(allspikes((allspikes>=(tstart-startoffset))& (allspikes<tstart)),...
                thisTrialSpikeTimes);
            thisTrialSpikeTimes = thisTrialSpikeTimes - tstart; % 0 = Trial ON of first template subtrial

            % modify spiketimes to account for the extra samples
            thisTrialSpikes = thisTrialSpikeTimes + startoffset; % just to start from zero
            NaNON = find(diff(isnan(Replay.TemplateTraces.Lever{whichreplay}))==1)+1;
            NaNOFF = find(diff(isnan(Replay.TemplateTraces.Lever{whichreplay}))==-1);
            if isnan(Replay.TemplateTraces.Lever{whichreplay}(1,1))
                NaNOFF(1,:) = [];
            end
            NaNPeriods = [ NaNON NaNOFF];
            % first flag out any spikes that would fall during these NaN
            % periods
            for k = 1:size(NaNPeriods)
                t1 = (NaNPeriods(k,1)-1)/SampleRate;
                t2 = NaNPeriods(k,2)/SampleRate;
                thisTrialSpikes((thisTrialSpikes>t1)&(thisTrialSpikes<=t2)) = NaN;
            end
            % now adjust spiketimes to account for extra samples
            for k = 1:size(NaNPeriods)
                t2 = NaNPeriods(k,2)/SampleRate;
                offset = (diff(NaNPeriods(k,:))+1)/SampleRate;
                thisTrialSpikes(thisTrialSpikes>t2) = thisTrialSpikes(thisTrialSpikes>t2) - offset;
            end

            [myPSTH,~,myRaster] = MakePSTH(thisTrialSpikes',0,...
                [0 1000*ceil(TemplateTrialLength)],'downsample',SampleRate,'kernelsize',smoothingfactor);
            PSTH(1,1:numel(myPSTH),i) = myPSTH;
            Raster(1,1:size(myRaster,2),i) = myRaster;

            if plotephysfigs || savereplayfigs
                % plot raster
                subplot(units_per_fig,2,Rasterplot);
                row_idx = 1;
                row_idx = (x*2)-1;
                PlotRaster(thisTrialSpikes,row_idx,'k');

%                 % plot FR
%                 subplot(units_per_fig,2,FRplot);
%                 plot((1/SampleRate)*(1:numel(myPSTH)),myPSTH,'k');
% 
%                 FRmax = max(FRmax,max(myPSTH));
            end

            % plot the replay spike times
            MyTrials = Replay.(TraceTag).TrialIDs{whichreplay};
            for thisTrial = 1:numel(MyTrials)

                % get odor valve ON-OFF times from OEPS
                OdorTTLs = Replay.TTLs.OdorValve{thisTrial}; % TS of odor valve ON-OFF w.r.t. Trial start TTL
                if size(OdorTTLs,1)>nsubtrials % some replays have an extra odor transition at the very beginning - to shut off odor from pre-replay trial
                    OdorTTLs(1,:) = [];
                end
                % Add one more column for TrialStart w.r.t. Odor ON
                OdorTTLs(:,end+1) = OdorTTLs(:,1) - TrialInfo.OdorStart(TemplateTrials,1);
                % force first subtrial to start at ~0
                OdorTTLs(1,end) = OdorTTLs(1,2) - FirstTrialDuration;
                % convert TS to real OEPS time base
                OdorTTLs(:,[1 2 5]) = OdorTTLs(:,[1 2 5]) + TTLs.Trial(MyTrials(thisTrial),1);

                tstart = OdorTTLs(1,5);
                t1 = tstart - startoffset;
                t2 = TTLs.Trial(MyTrials(thisTrial)+1,1); % until next trial start
                thisTrialSpikeTimes = allspikes((allspikes>=t1)& (allspikes<t2)) - t1;

                % FR: Use the original spiketimes to get PSTH, split the PSTH
                [myPSTH,~,myRaster] = MakePSTH(thisTrialSpikeTimes',0,...
                    +[0 1000*ceil(ReplayTrialLength)],...
                    'downsample',SampleRate,'kernelsize',smoothingfactor);

                PSTH(1+thisTrial,1:numel(myPSTH),i) = myPSTH';
                FRmax = max(FRmax,max(myPSTH));
                Raster(1+thisTrial,1:size(myRaster,2),i) = myRaster;

                if plotOL
                    if plotephysfigs || savereplayfigs
                        % plot raster
                        subplot(units_per_fig,2,Rasterplot);
                        row_idx = thisTrial+1;
                        row_idx = x*2;
                        PlotRaster(thisTrialSpikeTimes,row_idx,Plot_Colors('t')); %used to be r
                    end
                end
            end

        end

        %% Outputs
        %    Physiology = [];
        %     Physiology(whichreplay).PSTH = PSTH;
        %     Physiology(whichreplay).Correlation = PSTHCorr(PSTH,whichUnits);


    end

end