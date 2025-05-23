
function [MyTraces,timestamps,PSTH,Raster] = ProcessOpenLoopTrials(Replay, TrialInfo, SingleUnits, TTLs, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('whichreplays', [], @(x) isnumeric(x));
params.addParameter('plotfigures', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotephys', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('savefigures', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('savepdfs', false, @(x) islogical(x) || x==0 || x==1);
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
savereplaypdfs = params.Results.savepdfs;
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
global pdfPosition;
global MyPDFname;

if isempty(allreplays)
    allreplays = 1:numel(Replay.TemplateTraces.TrialIDs);
end

if ~isempty(SingleUnits) && isempty(whichUnits)
    whichUnits = 1:size(SingleUnits,2);
end

PassiveReplays = Replay.TTLs.TrialID(Replay.TTLs.TrialID>TrialInfo.TrialID(end));

for x = 1:numel(allreplays) % for every unique replay stretch
    whichreplay = allreplays(x);
    
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
    if plotreplayfigs || savereplayfigs
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
        
        if savereplayfigs
            figure(H1);
            saveas(gcf,[MyFileName,'_ReplayBehavior.fig']);
            close(gcf);
            
            figure(H2);
            saveas(gcf,[MyFileName,'_ReplayMotor.fig']);
            close(gcf);
        end
    end
    
    %% Ephys
    PSTH = []; % dimensions: [trial x samples x units]
    if ~isempty(SingleUnits)
        %units_per_fig = 5;
        if plotreplayfigs || savereplayfigs
            figure;
        end
        % get trial IDs for the closed loop template trials
        TemplateTrials = Replay.TemplateTraces.TrialIDs{whichreplay}; 
        nsubtrials = numel(TemplateTrials);
        ReplayTrialLength = size(Replay.ReplayTraces.Lever{whichreplay},1)/SampleRate; % in seconds
        TemplateTrialLength = size(Replay.TemplateTraces.Lever{whichreplay},1)/SampleRate; % in seconds
        FirstTrialDuration = TrialInfo.Duration(TemplateTrials(1));
        
        for i = 1:numel(whichUnits) % for every cell
            MyUnit = whichUnits(i);
            FRmax = 25; % for plotting
            
            if plotephysfigs || savereplayfigs
                if mod(i,units_per_fig)
                    FRplot = 2*mod(i,units_per_fig);
                else
                    FRplot = 2*units_per_fig;
                end
                Rasterplot = FRplot - 1;
                
                subplot(units_per_fig,2,Rasterplot);
                % plot the trial structure
                PlotBehavior(timestamps,[],[],[],[],Trial,[],numel(PassiveReplays)+trials_per_replay+1);
                set(gca,'YLim',...
                    [-0.4 numel(PassiveReplays)+trials_per_replay+1.4],...
                    'YTick',[0 5],'TickDir','out','XLim',[0 round(timestamps(end))], ...
                    'XTick', []);
                axis manual;
                hold on;
                
                subplot(units_per_fig,2,FRplot);
                % plot the trial structure
                PlotBehavior(timestamps,[],[],[],[],Trial,[],FRmax);
                set(gca,'YLim',[0 FRmax],'YTick',[0 50 100],'TickDir','out','XLim',[0 round(timestamps(end))], ...
                    'XTick', []);
                axis manual;
                hold on;
            end
            
            % Get Spikes
            allspikes = SingleUnits(MyUnit).spikes;
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
                PlotRaster(thisTrialSpikes,row_idx,'k');
                
                % plot FR
                subplot(units_per_fig,2,FRplot);
                plot((1/SampleRate)*(1:numel(myPSTH)),myPSTH,'k');
                
                FRmax = max(FRmax,max(myPSTH));
            end
            
            % plot the replay spike times
            MyTrials = Replay.ReplayTraces.TrialIDs{whichreplay};
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
                        PlotRaster(thisTrialSpikeTimes,row_idx,Plot_Colors('t')); %used to be r
                    end
                end
            end
            
            % plot avg. FR for the replay trials
            if plotOL
                if plotephysfigs || savereplayfigs
                    
                    subplot(units_per_fig,2,FRplot);
                    temp = squeeze(PSTH(2:(trials_per_replay+1),:,i));
                    meanFR = mean(temp);
                    semFR = std(temp)/sqrt(trials_per_replay);
                    if ShadeErrorbars
                        MyShadedErrorBar((1/SampleRate)*(1:size(PSTH,2)),meanFR,semFR,Plot_Colors('t'),[],0.5); %used to be r
                    else
                        plot((1/SampleRate)*(1:size(PSTH,2)),meanFR,'color',Plot_Colors('t'),'Linewidth',0.5);
                        plot((1/SampleRate)*(1:size(PSTH,2)),meanFR+semFR,'color',Plot_Colors('t'),'Linewidth',0.25);
                        plot((1/SampleRate)*(1:size(PSTH,2)),meanFR-semFR,'color',Plot_Colors('t'),'Linewidth',0.25);
                    end
                    %plot((1/SampleRate)*(1:size(PSTH,2)),mean(,1),'r');
                    
                    if isempty(PassiveReplays)
                        title(['Unit# ',num2str(MyUnit)]);
                        
                        if mod(i,units_per_fig) == 0
                            if savereplayfigs
                                saveas(gcf,[MyFileName,'_MyUnits_',num2str(MyUnit/units_per_fig),'.fig']);
                                close(gcf);
                            end
%                             if savereplaypdfs
%                                 set(gcf, 'Units', 'Normalized', 'OuterPosition', pdfPosition);
%                                 exportgraphics(gcf, ...
%                                     [MyFileName,'_MyUnits_',num2str(MyUnit/units_per_fig),'.pdf'],...
%                                     'ContentType','vector','Append',true);
%                                     close(gcf);
%                             end
                            figure;
                        end
                    end
                end
            end
            
            % plot the passive replays as well
            if ~isempty(PassiveReplays) && x==1
                nReplays = numel(MyTrials);
                MyTrials = PassiveReplays;
                for thisTrial = 1:numel(MyTrials)
                    
                    % get odor valve ON-OFF times from OEPS
                    OdorTTLs = Replay.TTLs.OdorValve{thisTrial+nReplays}; % TS of odor valve ON-OFF w.r.t. Trial start TTL
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
                    if MyTrials(thisTrial)>=size(TTLs.Trial,1)
                        t2 = TTLs.Trial(MyTrials(thisTrial),2) + startoffset;
                    else
                        t2 = TTLs.Trial(MyTrials(thisTrial)+1,1); % until next trial start
                    end
                    thisTrialSpikeTimes = allspikes((allspikes>=t1)& (allspikes<t2)) - t1;
                
                    [myPSTH,~,myRaster] = MakePSTH(thisTrialSpikeTimes',0,...
                        +[0 1000*ceil(ReplayTrialLength)],...
                        'downsample',SampleRate,'kernelsize',smoothingfactor);
                    
                    PSTH(1+trials_per_replay+thisTrial,1:numel(myPSTH),i) = myPSTH';
                    FRmax = max(FRmax,max(myPSTH));
                    Raster(1+trials_per_replay+thisTrial,1:size(myRaster,2),i) = myRaster;
                    
                    if plotPR
                        if plotephysfigs || savereplayfigs
                            % plot raster
                            subplot(units_per_fig,2,Rasterplot);
                            if plotOL
                                row_idx = thisTrial+1+trials_per_replay;
                            else
                                row_idx = thisTrial+1;
                            end
                            PlotRaster(thisTrialSpikeTimes,row_idx,Plot_Colors('r')); % used to be t
                        end
                    end
                    
                end
                
                % plot avg. FR for the passive replay trials
                if plotPR
                    if plotephysfigs || savereplayfigs
                        
                        subplot(units_per_fig,2,FRplot);
                        
                        temp = squeeze(PSTH((trials_per_replay+2):end,:,i));
                        meanFR = mean(temp);
                        semFR = std(temp)/sqrt(numel(PassiveReplays));
                        %MyShadedErrorBar((1/SampleRate)*(1:size(PSTH,2)),meanFR,semFR,Plot_Colors('r'),[],0.5); % used to be t
                        if ShadeErrorbars
                            MyShadedErrorBar((1/SampleRate)*(1:size(PSTH,2)),meanFR,semFR,Plot_Colors('r'),[],0.5); %used to be r
                        else
                            plot((1/SampleRate)*(1:size(PSTH,2)),meanFR,'color',Plot_Colors('r'),'Linewidth',0.5);
                            plot((1/SampleRate)*(1:size(PSTH,2)),meanFR+semFR,'color',Plot_Colors('r'),'Linewidth',0.25);
                            plot((1/SampleRate)*(1:size(PSTH,2)),meanFR-semFR,'color',Plot_Colors('r'),'Linewidth',0.25);
                        end
                        
                        %                     plot((1/SampleRate)*(1:size(PSTH,2)),mean(squeeze(PSTH((trials_per_replay+2):end,:,i)),1),...
                        %                         'Color',Plot_Colors('t'));
                        
                        % rescale plot if necessary
                        if max(get(gca,'YLim'))<FRmax
                            set(gca,'YLim', [0 5*ceil(FRmax/5)]);
                        end
                        
                        title(['Unit# ',num2str(MyUnit), '; Clust# ', num2str(SingleUnits(MyUnit).id),...
                            '; tet# ',num2str(SingleUnits(MyUnit).tetrode), '; Fp ', num2str(SingleUnits(MyUnit).ISIquality(1)), ', ISIFrac ',num2str(SingleUnits(MyUnit).ISIquality(2)) ]);
                        if mod(i,units_per_fig) == 0
                            if savereplayfigs
                                saveas(gcf,[MyFileName,'_MyUnits_',num2str(MyUnit/units_per_fig),'.fig']);
                                close(gcf);
                            end
                            if savereplaypdfs
                                set(gcf, 'Units', 'Normalized', 'OuterPosition', pdfPosition);
                                if i == units_per_fig
                                    exportgraphics(gcf, ...
                                    MyPDFname,...
                                    'ContentType','vector');
                                else
                                exportgraphics(gcf, ...
                                    MyPDFname,...
                                    'ContentType','vector','Append',true);
                                end
                                close(gcf);
                            end
                            figure;
                        end
                    end
                end
                
            end
            
            if savereplayfigs
                saveas(gcf,[MyFileName,'_MyUnits_',num2str(MyUnit/units_per_fig),'.fig']);
                close(gcf);
            end
%             if savereplaypdfs
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', pdfPosition);
%                 exportgraphics(gcf, ...
%                     [MyFileName,'_MyUnits_',num2str(MyUnit/units_per_fig),'.pdf'],...
%                     'ContentType','vector','Append',true);
%                 close(gcf);
%             end
            
        end
        
        %% Outputs
        %    Physiology = [];
        %     Physiology(whichreplay).PSTH = PSTH;
        %     Physiology(whichreplay).Correlation = PSTHCorr(PSTH,whichUnits);
        
        
    end
    
end