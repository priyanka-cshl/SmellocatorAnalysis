function [x, FR, BinOffset] = PlotTuningTrials(trialsdone, whichUnit, whichOdor, AlignedSpikes, TuningTTLs, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('selectlocation', [], @(x) isnumeric(x));
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;
whichlocations = params.Results.selectlocation;
psth = params.Results.psth;

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
whichUnit = abs(whichUnit);

thisUnitSpikes = AlignedSpikes(:,whichUnit);
whichodor = whichOdor;
% get the trial order
if ~isempty(whichlocations)
    whichTrials = intersect(find(TuningTTLs(:,5)==whichOdor),...
        find(round(TuningTTLs(:,7)/10)==whichlocations/10));
else
    whichTrials = find(TuningTTLs(:,5)==whichOdor);
end
whichTrials = [whichTrials TuningTTLs(whichTrials,7) TuningTTLs(whichTrials,8)]; %[tuningtrial# location #trialID]
whichTrials = sortrows(whichTrials,2); % sort by location

Xlims = [-1 2];
SpikesPSTH = [];

if plotting
    
    if plotevents
        % Plot Events
        %EventPlotter(myEvents);
        % Plot TrialType
        TrialTypePlotter(whichTrials(:,2),whichodor,Xlims,0);
    end
    
    % Plot Spikes
    for x = 1:size(whichTrials,1)
        if plotevents
            % Plot Odor On periods - adjust times if needed
            ZoneTimes = [0 (TuningTTLs(whichTrials(x,1),6) - TuningTTLs(whichTrials(x,1),4))];
            InZonePlotter(ZoneTimes', (x + trialsdone));
        end
        if plotspikes
            % Plot Spikes
            OdorStart = TuningTTLs(whichTrials(x,1),4);
            thisTrialSpikes = thisUnitSpikes.spikes(...
                find( (thisUnitSpikes.spikes>(OdorStart+Xlims(1))) & ...
                (thisUnitSpikes.spikes<(OdorStart+Xlims(2))) ) );
            thisTrialSpikes = thisTrialSpikes - OdorStart;
            PlotRaster(thisTrialSpikes,(x + trialsdone),Plot_Colors('r'));
            
            SpikesPSTH = vertcat(SpikesPSTH, [thisTrialSpikes (x)*ones(numel(thisTrialSpikes),1)]);
        end
    end
end

% calculate PSTH
FR = [];
BinOffset = round(Xlims(1)*1000);
if psth
    % aggregate psth of all trials
    if ~isempty(SpikesPSTH)
        FR(1,:) = VerySimplePSTH(SpikesPSTH(:,1), x, BinOffset,'downsample',500);
    end
end

% if psth
%     
%     
%     for i = -1:1:3
%         whichsniffs = find(AllSniffs(:,5)==i);
%         if ~isempty(whichsniffs)
%             %FR{i+2} = MakePSTH_v4(SpikesPlot(find(ismember(SpikesPlot(:,2),whichsniffs)),1),numel(whichsniffs),BinOffset,'downsample',500,'kernelsize',20);
%             FR{i+2} = MakeSniffTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),whichsniffs)),1),...
%                 AllSniffs(whichsniffs,15),...
%                 BinOffset,'downsample',500,'kernelsize',20);
%         else
%             FR{i+2} = [];
%         end
%     end
% else
%     FR = [];
% end

x = x + trialsdone;

% 
% for TZ = 1:12
%     thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1));
%     Events2Align = Offset(find(whichTrials(:,2)==TZ),1);
%     [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
%     AlignedFRs(TZ,1:numel(myFR)) = myFR;
%     RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
% end
% 
% x = size(allTrials,1);

%%
    function EventPlotter(myEvents)
        ticklength = 0.8;
        Y = ((1:size(myEvents,1))' + repmat([-ticklength 0 0],size(myEvents,1),1) )';
        
        % Odor ON
        X = [repmat(myEvents(:,1),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('r'),'Linewidth',2);
        
        % Trial OFF
        X = [repmat(myEvents(:,3),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('r'),'Linewidth',2);
        
        % Rewards
        X = [repmat(myEvents(:,2),1,2) NaN*ones(size(myEvents,1),1)]';
        plot(X(:),Y(:),'Color',Plot_Colors('b'),'Linewidth',2);
    end

    function InZonePlotter(TS, rowidx)
        foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.4);
        foo.EdgeColor = 'none';
        if ~isempty(TS)
            foo.Vertices = [ ...
                reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
                repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
            foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
        end
    end

    function TrialTypePlotter(TrialList,OdorType,Xlims,trialsdone)
        X = [Xlims(1) Xlims(2) Xlims(2) Xlims(1)];
        U = unique(TrialList);
        y1 = trialsdone;
        boxcolor(1,:) = Plot_Colors(['Odor',num2str(OdorType)]);
        boxcolor(2,:) = boxcolor-0.2;
        for j = 1:numel(U)
            y2 = y1 + numel(find(TrialList==U(j)));
            Y = [y1 y1 y2 y2];
            if trialsdone
                fill(X,Y,boxcolor(1,:),'EdgeColor','k');
            else
                fill(X,Y,boxcolor(1,:),'EdgeColor','none');
            end
            y1 = y2;
            boxcolor = flipud(boxcolor);
        end
    end
end