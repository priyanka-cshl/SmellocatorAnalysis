%% script to show sniff aligned cells - in ITI
% cells: Q3_20221019_r0_18, Q8_20221204_r0_49, S12_20230731_r0_10,
% S12_20230731_r0_58
warptype = 0;
sortbylocation = 0;
plotevents = 1;
plotspikes = 1;
psth = 1;
BinOffset = -100; % in ms
shade_exhalation = 1;

figure;
for examples = 1:4
    switch examples
        case 1
            SessionTag = 'S12_20230731_r0_processed.mat';   UnitID = 58;    OdorID = 3;
        case 2
            SessionTag = 'S12_20230731_r0_processed.mat';   UnitID = 31;    OdorID = 3;
        case 3
            SessionTag = 'S12_20230731_r0_processed.mat';    UnitID = 44;   OdorID = 3;
        case 4
            SessionTag = 'Q4_20221109_r0_processed.mat';    UnitID = 55;    OdorID = 2;        
    end

MouseName = regexprep(SessionTag,'_(\w+)_processed.mat','');
Paths = WhichComputer();
MySession = fullfile(Paths.ProcessedSessions,MouseName,SessionTag);

[TrialAligned, TrialInfo, ...
    ReplayAligned, ReplayInfo, ...
    TuningAligned, TuningInfo, ...
    AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

% get sniff time stamps and info for the sniffs we want to plot
[SelectedSniffs] = SelectSniffs_forKernelFits(TrialAligned, TrialInfo, [1 2 3]);

% some processing
SniffParams = [];
for whichodor = 1:3
    % get the spike raster
    [~, SniffParams_temp] = ...
        SniffAlignedFRs(SelectedSniffs{whichodor}, [], 'onlyparams', 1);
    SniffParams = vertcat(SniffParams, SniffParams_temp);
end
% sniff params
% 1-4: [currsniffstate currsniffloc currinhend currsniffend ...
% 5-6:  currsniffTrialID currsniffIndex ...
% 7:10: prevsniffstate prevsniffloc previnhstart previnhend
% 11:   actual sniff start time in OEPS timebase]
ITISniffs = SniffParams(find(SniffParams(:,1)==-1),:);
tokeep = randperm(size(ITISniffs,1),1000);
ITISniffs = ITISniffs(tokeep,:);
% sort by sniff duration
ITISniffs = sortrows(ITISniffs,4);

OdorSniffs = SniffParams(find(SniffParams(:,1)==OdorID),:);
tokeep = randperm(size(OdorSniffs,1),1000);
OdorSniffs = OdorSniffs(tokeep,:);
% sort by sniff duration
OdorSniffs = sortrows(OdorSniffs,4);

% switch examples
%     case 1 % 4277 sniffs
%         % remove the last 177
%         tokeep = sort(randperm(4100,1000));
%     case 2 % 4277 sniffs
%         tokeep = sort( [(3000 + randperm(1277,1000)) randperm(3000,1000)]);
%     case 3 % 2207 sniffs
%         % remove the shortest 107 and longest 100
%         tokeep = sort( 25 + randperm(2180,2000));
%         %tokeep = 108:2107;
%     case 4 % 3019 sniffs
%         tokeep = sort(randperm(3019,2000));
% end
%ITISniffs = ITISniffs(tokeep,:);

if sortbylocation
    % sort by previous sniff duration
    OdorSniffs = sortrows(OdorSniffs,2,'descend'); % previous sniff durations are negative
end

sniffmax(examples)  = max(ITISniffs(:,4)) + 0.2;
thisUnitSpikes = AllUnits.Spikes{UnitID};

AllSniffs = [ITISniffs; OdorSniffs];

%% Plotting sniff periods and spikes
subplot(1,4,examples)
hold on;

SpikesPlot = [];
for x = 1:size(AllSniffs,1)    
    if plotevents
        % Plot inhalation periods - adjust times if needed
        InhalationTimes = [AllSniffs(x,[9 10]) 0 AllSniffs(x,3)];
        InhalationTimes = reshape(InhalationTimes,2,2)';
        
        ExhalationTimes = [AllSniffs(x,10) 0 AllSniffs(x,3:4)];
        ExhalationTimes = reshape(ExhalationTimes,2,2)';
        
        if warptype % warp by inhalation duration
            InhalationTimes = InhalationTimes * (mean(AllSniffs(:,3))/AllSniffs(x,3));
            ExhalationTimes = ExhalationTimes * (mean(AllSniffs(:,3))/AllSniffs(x,3));
        end
        if ~shade_exhalation
            SniffPlotter(InhalationTimes', x, Plot_Colors('TZ'));
        else
            if AllSniffs(x,1) < 0
                SniffPlotter(ExhalationTimes', x, Plot_Colors('Odor3'));
            else
                SniffPlotter(ExhalationTimes', x, Plot_Colors('Odor2'));
            end
                
        end
    end
    
    if plotspikes || psth
        ts = [AllSniffs(x,9) 0 AllSniffs(x,4) AllSniffs(x,3)] + AllSniffs(x,11); 
        % [prevsniffstart thisniffstart thissniffend thisinhend]        
        whichSpikes = intersect(find(thisUnitSpikes>=ts(1)), find(thisUnitSpikes<=ts(3)));
        
        % align to inhalation start
        thisSniffSpikes = thisUnitSpikes(whichSpikes) - ts(2);
        tmax = ts(3) - ts(2);
        
        if warptype
            thisSniffSpikes = thisSniffSpikes * (mean(AllSniffs(:,3))/AllSniffs(x,3));
            tmax = tmax * (mean(AllSniffs(:,3))/AllSniffs(x,3));
        end
        
        SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
        
    end
end
nSniffs(examples) = size(AllSniffs,1);

% plot spikes
plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize', 2);
%set(gca,'YLim', [0 2200], 'XLim',[-0.05 0.75],'TickDir','out');
set(gca,'XLim',[-0.05 0.75],'TickDir','out');

% PSTH calculation
if psth
    FR{examples,1} = MakeSniffTriggeredPSTH(SpikesPlot(:,1),AllSniffs(:,4),...
                BinOffset,'downsample',500,'kernelsize',20);
    % also make psths for the top half and bottom half?    
    whichsniffs = 1:1000;
    FR{examples,2} = MakeSniffTriggeredPSTH(...
                        SpikesPlot(find(ismember(SpikesPlot(:,2),whichsniffs)),1),...
                        AllSniffs(whichsniffs,4),...
                        BinOffset,'downsample',500,'kernelsize',20);
    whichsniffs = 1000 + (1:1000);
    FR{examples,3} = MakeSniffTriggeredPSTH(...
                        SpikesPlot(find(ismember(SpikesPlot(:,2),whichsniffs)),1),...
                        AllSniffs(whichsniffs,4),...
                        BinOffset,'downsample',500,'kernelsize',20);  
end

end

set(gcf, 'Position', [281 405 1160 392]);
set(gcf,'Renderer','painters');
if ~sortbylocation
    print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/OdorsniffsExampleRasters.eps',...
    '-depsc','-tiff','-r300','-painters');
else
    print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/OdorsniffsExampleRasters_LocationSorted.eps',...
    '-depsc','-tiff','-r300','-painters');
end


%% plot the psths
figure;
for examples = 1:4
    subplot(1,4,examples);
    ts = (1:length(FR{examples,1}))*0.002;
    ts = ts + (BinOffset/1000);
    plot(ts,FR{examples,1},'k');
    hold on
    ts = (1:length(FR{examples,2}))*0.002;
    ts = ts + (BinOffset/1000);
    plot(ts,FR{examples,2},'b');
    ts = (1:length(FR{examples,3}))*0.002;
    ts = ts + (BinOffset/1000);
    plot(ts,FR{examples,3},'r');
    set(gca,'XLim',[-0.05 0.75],'YLim',[0 30],'TickDir','out');   
end

set(gcf, 'Position', [281 405 1160 150]);
set(gcf,'Renderer','painters');
if ~sortbylocation
    print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/OdorsniffsExampleRasters_FRs.eps',...
        '-depsc','-tiff','-r300','-painters');
else
    print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/OdorsniffsExampleRasters_LocationSorted_FRs.eps',...
        '-depsc','-tiff','-r300','-painters');
end

%%
%handles.(['spikesplot',num2str(i)]) = plot(NaN, NaN, '.k','Markersize', 0.5);
function SniffPlotter(AllTS, rowidx, boxcolor)
%foo = fill(NaN,NaN,Plot_Colors('TZ'),'FaceAlpha',0.8);
for t = 1:size(AllTS,2)
    TS = AllTS(:,t);
    foo = fill(NaN,NaN,boxcolor,'FaceAlpha',0.7);
    foo.EdgeColor = 'none';
    if ~isempty(TS)
        foo.Vertices = [ ...
            reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
            repmat([rowidx-1 rowidx rowidx rowidx-1]',size(TS,2),1)];
        foo.Faces = reshape(1:2*numel(TS),4,size(TS,2))';
    end
end
end
