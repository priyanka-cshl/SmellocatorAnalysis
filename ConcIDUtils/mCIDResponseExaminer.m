function [] = mCIDResponseExaminer() %(myDir,myStimFile)

myKsDir = '/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31/';
myStimFile = '/mnt/data/Sorted/T3/250516/250516_13_55.txt';

% get trial and odor valve transitions
temp = load(fullfile(myKsDir,'myTTLfile_1.mat'));
ch = temp.TTLs.data;
states = temp.TTLs.info.eventId;
EventTS = temp.TTLs.timestamps;
EventTS = EventTS - temp.TTLs.offset; % subtract start timestamp to align with spikes

% Trial timestamps
TTLs.Trial(:,1) = EventTS(ch==1 & states,1);
while size(EventTS(ch==1 & ~states,1),1) > size(TTLs.Trial,1)
    TTLs.Trial(end,:) = [];
end
TTLs.Trial(:,2) = EventTS(ch==1 & ~states,1);
TTLs.Trial(1,:) = [];

% Odor timestamps
TTLs.Odor(:,1) = EventTS(ch==2 & states,1);
while size(EventTS(ch==2 & ~states,1),1) > size(TTLs.Odor,1)
    TTLs.Odor(end,:) = [];
end
TTLs.Odor(:,2) = EventTS(ch==2 & ~states,1);
TTLs.Odor(1,:) = [];

% Get stimulus info from the odor machine file
StimFile = readmatrix(myStimFile);
StimSettings.timing = StimFile(1:6)'; % [whichmachine(0=16 odors) pre-stim stim post-stim iti reps]
StimFile(1:6,:) = [];

% some info about number of stimuli and various reps etc
nStim = unique(StimFile);
nTypes = 2; % no. of concentrations checked
nReps = StimSettings.timing(end);
extraStim = unique(StimFile(((numel(nStim)*nReps)+1):end));
extraReps = (length(StimFile) - (numel(nStim)*nReps))/numel(extraStim);
nStim = nStim(~ismember(nStim,extraStim));
StimSettings.miniOdors = nStim;
StimSettings.miniReps = nReps;
StimSettings.megaOdors = extraStim;
StimSettings.megaReps = extraReps;

if size(TTLs.Trial,1) == size(StimFile,1)
    TTLs.Trial(:,3) = TTLs.Trial(:,2) - TTLs.Trial(:,1); % Trial duration
    % add odor info to Trial TTLs
    TTLs.Trial(:,4) = StimFile; % odor identity
else
    keyboard;
end

% Add odor starts and stops
for t = 1:size(TTLs.Trial,1)
    whichrows = find(TTLs.Odor(:,1)>TTLs.Trial(t,1) & TTLs.Odor(:,2)<TTLs.Trial(t,2));
    odorTTLs = TTLs.Odor(whichrows,:)';
    TTLs.Trial(t,6+(1:numel(odorTTLs))) = odorTTLs(:);
end

if ~exist(fullfile(myKsDir,'quickprocessOdorTTLs.mat'))
    save(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings');
end

%% Load spikes
SingleUnits = GetSingleUnits(myKsDir, 3);
nUnits = size(SingleUnits,2);

%% figures related
savefigs = 1;
mycolors = brewermap(10,'YlOrRd');
mycolors(1:2,:) = [];
nRows = 4; nCols = 6;
conc = 1;

for n = 1:nUnits
    FigureName = ['Unit ',num2str(n)]; % one figure per cell
    figure('Name',FigureName);
    thisUnitSpikes = SingleUnits(n).spikes;
    for odor = 1:numel(nStim)
        if odor > 5
            subplot(nRows,nCols,odor+1);
        else
            subplot(nRows,nCols,odor);
        end
        hold on

        whichTrials = find(TTLs.Trial(:,4)==nStim(odor));
        odorON = []; SpikesPlot = [];
        for rep = 1:numel(whichTrials)
            ts = TTLs.Trial(whichTrials(rep),[1 2 7 8 9 10]); % trial start, stop, odor start, stop, purge start, stop
            thistrialspikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<ts(2)));
            thistrialspikes = thisUnitSpikes(thistrialspikes) - ts(3);
            odorON = vertcat(odorON, [ts(3:6)-ts(3) , rep]);
            SpikesPlot = vertcat(SpikesPlot, [thistrialspikes rep*ones(numel(thistrialspikes),1) ]);
        end
        
        plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 3, 'color', 'k'); %mycolors(conc*2,:));
        set(gca,'XTick',[],'YTick',[],'XLim',[-5 20], 'YLim',[0 (nReps+1)/500]);
        for o = 1:4
            plot(odorON(:,o),odorON(:,end)/500,'r');
        end
    end
    
    for odor = 1:numel(extraStim)
        whichplots = 12 + (odor:6:12);
        subplot(nRows,nCols,whichplots);
        hold on

        whichTrials = find(TTLs.Trial(:,4)==extraStim(odor));
        odorON = []; SpikesPlot = [];
        for rep = 1:numel(whichTrials)
            ts = TTLs.Trial(whichTrials(rep),[1 2 7 8 9 10]); % trial start, stop, odor start, stop, purge start, stop
            thistrialspikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<ts(2)));
            thistrialspikes = thisUnitSpikes(thistrialspikes) - ts(3);
            odorON = vertcat(odorON, [ts(3:6)-ts(3) , rep]);
            SpikesPlot = vertcat(SpikesPlot, [thistrialspikes rep*ones(numel(thistrialspikes),1) ]);
        end
        
        plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 3, 'color', 'k'); %mycolors(conc*2,:));
        set(gca,'XTick',[],'YTick',[],'XLim',[-5 20], 'YLim',[0 (extraReps+1)/500]);
        for o = 1:4
            plot(odorON(:,o),odorON(:,end)/500,'r');
        end
    end

    set(gcf,'Position',[813 752 1108 217]);
    if savefigs
        set(gcf,'Color','w');
        set(gcf,'renderer','Painters');
        figPosition = [0.2089    0.6907    0.2911    0.2796];
        set(gcf, 'Units', 'Normalized', 'OuterPosition', figPosition);
        if n == 1
            exportgraphics(gcf, ...
                fullfile(myKsDir,'OdorMaps','OdorSummary.pdf'),...
                'ContentType','vector');
        else
            exportgraphics(gcf, ...
                fullfile(myKsDir,'OdorMaps','OdorSummary.pdf'),...
                'ContentType','vector','Append',true);
        end
        close(gcf);
    end
end

end