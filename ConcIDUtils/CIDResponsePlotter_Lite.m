function [] = CIDResponsePlotter_Lite(sessionPath,mySpikes,unitID)
if nargin < 3
    unitID = nan;
end

StimSettings.SessionType = '16_Odors';
%sessionPath = '/mnt/data/CID/forGLM/E3_2022-07-07_16-47-23_cid_processed.mat'; % path to the preprocessed session 
load(sessionPath);
% make TTLs
OdorValve = TracesOut.Odor{1};
OdorValve(OdorValve<0) = 0;

OdorTTls = find(diff(OdorValve)>0) + 1; % On
OdorTTls(:,2) = find(diff(OdorValve)<0);
OdorIDs = OdorValve(OdorTTls(:,1));
OdorTTls = TracesOut.Timestamps{1}(OdorTTls);
odorDuration = round(median(OdorTTls(:,2)-OdorTTls(:,1)))*1000;

nStim = unique(OdorIDs);
align2sniffs = 0;

%% Actual Plotting
if strcmp(StimSettings.SessionType,'newCID')
    if bothPulses
        plotWidth = [-4 (odorDuration*5)/1000];
    else
        plotWidth = [-6 (odorDuration*2)/1000];
    end
else
    plotWidth = (odorDuration*2)/1000*[-0.5 1];
end

switch StimSettings.SessionType
    case {'newCID', '16_Odors'}
        mycolors = brewermap(numel(nStim),'Dark2');
        foo = round(numel(nStim)/2)+1;
        mycolors(foo:end,:) = mycolors(foo:end,:)/2; 

    case {'16_Concs*'}
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        for c = 0:0.2:0.6
            lighter = basecolors + (1 - basecolors)*c;
            mycolors = vertcat(lighter, mycolors);
        end
        StimSettings.SessionType = '16_Odors';

    case '16_Concs'
        basecolors = brewermap(5,'Dark2');
        mycolors = [];
        scale = 0.6:-0.2:0;
        for c = 1:4
            mycolors(:,:,c) = basecolors + (1 - basecolors)*scale(c);
        end
        mycolors = reshape(permute(mycolors, [3 1 2]), 20, 3);      % rows: [color1's N shades; color2's N shades; ...]
        StimSettings.SessionType = '16_Odors';

    case 'otherwise'
        nConcs = numel(nTypes);
        extraColors = 2;
        mycolors = brewermap(2*(nConcs+extraColors),'YlOrRd');
        mycolors(1:(2*extraColors),:) = []; % too light
        trialsPerOdor = size(TTLs.Trial,1)/nStim;
        unitsPerPlot = 5;
        plots_rows = unitsPerPlot;
        plots_cols = numel(nStim);
        panelsPerPlot = plots_rows*plots_cols; % one per unit
        %plotWidth = [-6 (StimSettings.timing(3)*2)/1000];
        plotWidth = (odorDuration*2)/1000*[-0.5 1];
end

%% Actual Plotting - input spikes
        figure;
for foo = 1:(1+~isnan(unitID))
    if ~isnan(unitID)
        subplot(1,2,foo);
        hold on;
        if foo == 2
            whichUnit = find([SingleUnits.id]==unitID);
            mySpikes = SingleUnits(whichUnit).spikes;
        end
    else

        hold on;
    end
    % if 16 odors, loop by odors per unit for one subplot,
    odorON = [];
    repsDone = 0;
    for odor = 1:numel(nStim)
        SpikesPlot = [];
        whichtrials = find(OdorIDs==nStim(odor));
        nreps = numel(whichtrials);
        for rep = 1:numel(whichtrials)
            thisTrial = whichtrials(rep);
            t = OdorTTls(thisTrial,:) + [-odorDuration odorDuration]/1000;
            thistrialspikes = [];
            thistrialspikes = mySpikes(mySpikes>=t(1) & mySpikes<=t(2)) - OdorTTls(thisTrial,1);
            SpikesPlot = vertcat(SpikesPlot, ...
                [thistrialspikes ...
                (rep + repsDone)*ones(numel(thistrialspikes),1)]...
                );
            odorON = vertcat(odorON, [0 , (rep + repsDone)]);
        end
        repsDone = repsDone + nreps;
        plot(SpikesPlot(:,1), SpikesPlot(:,2)/500, '.k','Markersize', 2, 'color', mycolors(odor,:));
    end

    plotHeight = (max(odorON(:,2))+1)/500;
    set(gca,'XTick',[],'YTick',[],'XLim',plotWidth,'YLim',[0 plotHeight]);
    plot(odorON(:,1),odorON(:,2)/500,'k');
    plot(odorON(:,1)+odorDuration/1000,odorON(:,2)/500,'k');
end

end
