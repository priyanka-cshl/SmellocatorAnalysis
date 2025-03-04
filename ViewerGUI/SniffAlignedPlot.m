function [nSniffs,FR,SpikesPlot] = SniffAlignedPlot(AllSniffs, thisUnitSpikes, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('alignto', 1, @(x) isnumeric(x)); % 1 - inhalation start, 2 - inhalation end
params.addParameter('warptype', 0, @(x) isnumeric(x)); % 0 - no warp, 1 - by sniff duration, 2 - by inhalation duration
params.addParameter('psthoffset', -1000, @(x) isnumeric(x));
params.addParameter('spikeplothandle', 1, @(x) ishghandle(x));

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;
psth = params.Results.psth;
alignto = params.Results.alignto;
warptype = params.Results.warptype;
BinOffset = params.Results.psthoffset; % timewindow before t0 included in the PSTH 
spikehandle = params.Results.spikeplothandle;

% plotting sniff color
if plotevents
    SniffColors = colormap(brewermap(220,'Spectral'));
end
SpikesPlot = []; SpikesPSTH = [];

% Plot Spikes
for x = 1:size(AllSniffs,1)
    
    if plotevents
        % Plot exhalation periods - adjust times if needed
        ExhalationTimes = AllSniffs(x,[6 7 8 9 10 11]) - AllSniffs(x,alignto+6);
        ExhalationTimes = reshape(ExhalationTimes,2,3)';
        
        if warptype
            ExhalationTimes = ExhalationTimes * (mean(AllSniffs(:,14+warptype))/AllSniffs(x,14+warptype));
        end
        
        thissnifflocation = floor(AllSniffs(x,[12 13 14])+110);
        if any(isinf(thissnifflocation))
            fprintf('warning: Infs in sniff location\n');
            thissnifflocation(find(isinf(thissnifflocation))) = 1;
        end
        SniffPlotter(ExhalationTimes', x, SniffColors(thissnifflocation,:));
        
        if x == size(AllSniffs,1)
            line([-0.5 1], [x x], 'Color', 'k');
        end
    end
    
    if plotspikes || psth
        
        ts = AllSniffs(x,[5 7 11 9 8]); % [prevsniffstart thisniffstart nextsniffend thissniffend thisinhend]
        whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(3)));
        
        if alignto == 1 % align to inhalation start
            thisSniffSpikes = thisUnitSpikes(whichSpikes) - ts(2);
            tmax = ts(4) - ts(2);
        end
        if alignto == 2 % align to inhalation end
            thisSniffSpikes = thisUnitSpikes(whichSpikes) - ts(5);
            tmax = ts(4) - ts(5);
        end
        
        if warptype
            thisSniffSpikes = thisSniffSpikes * (mean(AllSniffs(:,14+warptype))/AllSniffs(x,14+warptype));
            tmax = tmax * (mean(AllSniffs(:,14+warptype))/AllSniffs(x,14+warptype));
        end
        
        SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
        
        if psth
            % for psth - ignore spikes that happened after the current sniff end
            thisSniffSpikes(thisSniffSpikes>tmax,:) = [];
            SpikesPSTH = vertcat(SpikesPSTH, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
        end
    end
end
nSniffs = size(AllSniffs,1);

% % Spike plotting
% if plotspikes && ~isempty(SpikesPlot)
%     %plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize',0.5);
%     set(spikehandle,'XData',SpikesPlot(:,1),'YData',SpikesPlot(:,2));
% end

% PSTH calculation
if psth
    % calculate by sniff type 
    for snifftype = -1:1:4
        whichsniffs = find(AllSniffs(:,4)==snifftype);
        if ~isempty(whichsniffs)
            SniffDurations = AllSniffs(whichsniffs,15);
            if warptype
                SniffDurations = (SniffDurations./AllSniffs(whichsniffs,14+warptype)) * ...
                 mean(AllSniffs(:,14+warptype));
            end
            FR{snifftype+2} = MakeSniffTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),whichsniffs)),1),...
                SniffDurations,...
                BinOffset,'downsample',500,'kernelsize',20);
        else
            FR{snifftype+2} = [];
        end
    end
else
    FR = [];
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

    function SniffPlotter(AllTS, rowidx, boxcolor)
        %foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.8);
        for t = 1:size(AllTS,2)
            TS = AllTS(:,t);
            foo = fill(NaN,NaN,boxcolor(t,:),'FaceAlpha',0.8);
            foo.EdgeColor = 'none';
            if ~isempty(TS)
                foo.Vertices = [ ...
                    reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
                    repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
                foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
            end
        end
    end
end