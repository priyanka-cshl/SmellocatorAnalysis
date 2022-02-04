MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up

%% get the data loaded
[TracesOut, ColNames, TrialInfo, SingleUnits, TTLs, SampleRate, TimestampAdjuster] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

%% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,size(TrialInfo.TrialID,2));

%% plot a figure, with 3 cols, one for each odor, sort trials by target zone type, and further by trial duration

for whichUnit = 58 %1:size(AlignedSpikes,2)
    thisUnitSpikes = AlignedSpikes(:,whichUnit);
    for i = 1:3 % 3 odors
        
        % get the trial sorting order
        whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
            find(TrialInfo.Odor==i));
        whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) TrialInfo.Duration(whichTrials)]; %#ok<AGROW>
        whichTrials = sortrows(whichTrials,2);
        for tz = 1:12
            whichTrials(whichTrials(:,2)==tz,:) = sortrows(whichTrials(whichTrials(:,2)==tz,:),3);
        end
        
        subplot(1,3,i); hold on
        
        % Plot TrialType
        TrialTypePlotter(whichTrials(:,2),i);
        
        % Plot all events
        EventPlotter(Events(whichTrials(:,1),:));
                
        % Plot Spikes
        for x = 1:size(whichTrials,1)
            
            % Plot Target Zone periods
            InZonePlotter(TrialInfo.InZone{whichTrials(x)}', x);
            
            % Plot Spikes
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
            PlotRaster(thisTrialSpikes,x,Plot_Colors('k'));
        end
        
        set(gca, 'XLim', [-1.2 6]);
        
        trialcounts(i) = x; %#ok<SAGROW>
    end
end

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

function TrialTypePlotter(TrialList,OdorType)
X = [-1.2 -1 -1 -1.2];
U = unique(TrialList);
    y1 = 0;
    boxcolor(1,:) = Plot_Colors(['Odor',num2str(OdorType)]);
    boxcolor(2,:) = boxcolor-0.2;
    for x = 1:numel(U)
        y2 = y1 + numel(find(TrialList==U(x)));
        Y = [y1 y1 y2 y2]; 
        fill(X,Y,boxcolor(1,:),'EdgeColor','none');
        y1 = y2;
        boxcolor = flipud(boxcolor);
    end
end
