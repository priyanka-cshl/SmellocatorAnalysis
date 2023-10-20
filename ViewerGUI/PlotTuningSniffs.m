function [x,FR,BinOffset] = PlotTuningSniffs(whichUnit, whichOdor, Spikes, TuningSniffs, TuningInfo, Sniffsdone, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotspikes', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('plotevents', true, @(x) islogical(x) || x==0 || x==1);
params.addParameter('psth', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('sortorder', 0, @(x) isnumeric(x)); % 0 - sniff duration, 1 - inhalation duration
params.addParameter('alignto', 1, @(x) isnumeric(x)); % 1 - inhalation start, 2 - inhalation end
params.addParameter('warptype', 0, @(x) isnumeric(x)); % 0 - no warp, 1 - by sniff duration, 2 - by inhalation duration
params.addParameter('selectlocation', [], @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
plotspikes = params.Results.plotspikes;
plotevents = params.Results.plotevents;
psth = params.Results.psth;
sortby = params.Results.sortorder;
alignto = params.Results.alignto;
warptype = params.Results.warptype;
whichlocations = params.Results.selectlocation;

plotting = whichUnit>0; % hack to use the same function for UnitViewer and for analysis
whichUnit = abs(whichUnit);

Xlims = -[0.1 1.1];

thisUnitSpikes = Spikes(whichUnit).spikes';
whichodor = whichOdor;

% which trials to use
whichtrials = find(TuningInfo(:,2)==whichOdor+1);
% get the sniffs for trials of that particular odor                        
whichsniffs = find(ismember(TuningSniffs(:,17),whichtrials));

% assemble a list of sniffs 
% [inh-start inh-end next-inh sniffID snifftype sniff-duration inh-duration Trial ID]

AllSniffs = TuningSniffs(whichsniffs,:);

if ~isempty(whichlocations)
    % also pull out tuning sniffs that share the same location as halts
    halt_sniffs = intersect(find((AllSniffs(:,5)>0)&(AllSniffs(:,5)<4)), ...
                    find(round(AllSniffs(:,6)/10)==whichlocations/10));
                
    AllSniffs = AllSniffs(halt_sniffs,:);
    
end

switch sortby
    case 0
        % sort sniff List by Sniff Type, then sniff duration, then inh duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[5 15 16 17]);
    case 1
        % sort sniff List by Sniff Type, then inh duration, then sniff duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[5 16 15 17]);
    case 2
        AllSniffs(find(AllSniffs(:,5)==0),5) = 1; 
        AllSniffs(find(AllSniffs(:,5)==1),5) = 1;
        AllSniffs(find(AllSniffs(:,5)==2),5) = 1; 
        % sort by sniff type, then odor location (col 6), then sniff duration, then inh duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[5 6 16 15 17]);
end

% plotting sniff color
SniffColors = colormap(brewermap(220,'Spectral'));

AllSniffs(find(isnan(AllSniffs(:,5))),:) = [];

SpikesPlot = []; SpikesPSTH = [];
if plotting
    
    % Plot Spikes
    for x = 1:size(AllSniffs,1)
        if plotevents
            % Plot Target Zone periods - adjust times if needed
            ExhalationTimes = AllSniffs(x,[8 9 2 3 12 13]) - AllSniffs(x,alignto);
            ExhalationTimes = reshape(ExhalationTimes,2,3)';
            if warptype
                ExhalationTimes = ExhalationTimes * (mean(AllSniffs(:,14+warptype))/AllSniffs(x,14+warptype));
            end
            thissnifflocation = floor(AllSniffs(x,6)+110);
            SniffPlotter(ExhalationTimes', x + Sniffsdone, SniffColors(thissnifflocation,:));
        end
        
        if plotspikes || psth
              
              s1 = find(thisUnitSpikes>=(AllSniffs(x,1)-2),1,'first'); % 2 seconds before the current sniff
              s2 = find(thisUnitSpikes<=(AllSniffs(x,3)+2),1,'last'); % 2 seconds after the current sniff
              
              switch alignto
                  case 1 % inhalation start
                      thissniffspikes = thisUnitSpikes(1,s1:s2) - AllSniffs(x,1);
                  case 2 % inhalation end
                      thissniffspikes = thisUnitSpikes(1,s1:s2) - AllSniffs(x,2);
              end
              
              if warptype
                  thissniffspikes = thissniffspikes * (mean(AllSniffs(:,14+warptype))/AllSniffs(x,14+warptype));
              end
              
              SpikesPlot = vertcat(SpikesPlot, [thissniffspikes' (x + Sniffsdone)*ones(numel(thissniffspikes),1)]);

              % for plotting PSTH
              switch alignto
                  case 1 % inhalation start
                      thissniffspikes(:,thissniffspikes>AllSniffs(x,15)) = [];
                  case 2 % inhalation end
                      thissniffspikes(:,thissniffspikes>(AllSniffs(x,15)-AllSniffs(x,16))) = [];
              end
              SpikesPSTH = vertcat(SpikesPSTH, [thissniffspikes' (x)*ones(numel(thissniffspikes),1)]);

        end
    end
end

BinOffset = -1000;
if plotspikes & ~isempty(SpikesPlot)
    plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize',0.5);
end

x = size(AllSniffs,1);
if psth
    for i = -1:1:3
        whichsniffs = find(AllSniffs(:,5)==i);
        if ~isempty(whichsniffs)
            %FR{i+2} = MakePSTH_v4(SpikesPlot(find(ismember(SpikesPlot(:,2),whichsniffs)),1),numel(whichsniffs),BinOffset,'downsample',500,'kernelsize',20);
            FR{i+2} = MakeSniffTriggeredPSTH(SpikesPSTH(find(ismember(SpikesPSTH(:,2),whichsniffs)),1),...
                AllSniffs(whichsniffs,15),...
                BinOffset,'downsample',500,'kernelsize',20);
        else
            FR{i+2} = [];
        end
    end
else
    FR = [];
end

x = x + Sniffsdone;

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

    function SniffPlotter(TS, rowidx, boxcolor)
        %foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.8);
        foo = fill(NaN,NaN,boxcolor,'FaceAlpha',0.8);
        foo.EdgeColor = 'none';
        if ~isempty(TS)
            foo.Vertices = [ ...
                reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
                repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
            foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
        end
    end

    function SniffTypePlotter(TrialList,OdorType,Xlims,trialsdone)
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