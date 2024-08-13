%% script to analyze the stereotypy in sniff waveforms
savefigs = 1;

SessionName = 'S12_20230731_r0';
% load the processed output from Smellocator
MySession = fullfile('/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/', ...
    SessionName(1:regexp(SessionName,'_','once')-1), ...
    [SessionName,'_processed.mat']);

if savefigs
    FigPath = ['/home/priyanka/Desktop/sniffPSTHPredictions/', SessionName(1:regexp(SessionName,'_','once')-1)];
    FigPath = fullfile(FigPath,'SniffInfluence3');
    if ~exist(FigPath,'dir')
        mkdir(FigPath);
        fileattrib(FigPath, '+w','a');
    end
end
% reprocess sniff traces on a trialwise basis
load(MySession,'Traces','TrialInfo','SniffTS');
Traces.OdorLocation     = Traces.Motor;
[SniffTimeStamps] = TrialWiseSniffs(TrialInfo,Traces); % [sniffstart sniffstop nextsniff odorlocation sniffslope stimstate trialID]
% remove overlapping sniffs
Sniffs = SniffTimeStamps(find(SniffTimeStamps(:,end)>0),:);

%% for converting snifftimestamps to OEPS base
load(MySession,'TimestampAdjust', 'TTLs');
% check that there was no clock drift
if any(abs(TTLs.Trial(1:numel(TrialInfo.Odor),2) - (TrialInfo.SessionTimestamps(:,2) + TimestampAdjust.ClosedLoop))>0.04)
    disp('clock drift in ephys and behavior files');
    keyboard;
else
    SniffsOEPS = Sniffs;
    % use TimestampAdjust to convert to OEPS timebase
    SniffsOEPS(:,1:3) = SniffsOEPS(:,1:3) + TimestampAdjust.ClosedLoop;
    % add a column for sniff duration
%    SniffsOEPS(:,7) = SniffsOEPS(:,3)-SniffsOEPS(:,1);
    SniffsOEPS(:,8) = SniffsOEPS(:,2)-SniffsOEPS(:,1);
    
    % make a continuous Odor Location Trace to extract odor location info in each sniff aligned bin
    load(MySession,'Traces');
    [TracesOut, whichtraces] = ConcatenateTraces2Mat(Traces);
    TraceTimestamps = TracesOut(:,find(strcmp(whichtraces,'Timestamps'))) + TimestampAdjust.ClosedLoop;
    OdorLocations   = TracesOut(:,find(strcmp(whichtraces,'Motor')));

    % also get a continuous sniff trace
    [SniffTrace] = FilterThermistor(TracesOut(:,find(strcmp(whichtraces,'Sniffs'))));

end

%% make a continuous Odor Location Trace to extract odor location info in each sniff aligned bin

%% analyzing spike rate vs sniff variability
load(MySession,'SingleUnits','Traces');
% SingleUnits(whichunit).spikes; % raw timestamps in OEPS base

binsize = 0.05; % in seconds
windowsize = 0.200; % in seconds
BinLims = [];
BinLims(:,1) = 0:binsize:windowsize;
BinLims(:,2) = circshift(BinLims(:,1),-1);
BinLims(end,:) = [];

ChosenUnits = [6 18 25 29 39 44 58 79 85];

for WhichUnit = 1:numel(ChosenUnits)

    SelectedUnit = ChosenUnits(WhichUnit); %39;

    BinnedPSTH   = [];
    for s = 1:size(SniffsOEPS,1) % every sniff
        thisSniffBins = BinLims + SniffsOEPS(s,1);
        for n = 1:size(BinLims,1)
            thisBinSpikes = find(SingleUnits(SelectedUnit).spikes>=thisSniffBins(n,1) & SingleUnits(SelectedUnit).spikes<thisSniffBins(n,2));
            BinnedPSTH(s,n) = numel(thisBinSpikes);

            % get the motor location in that bin
            x1 = find(TraceTimestamps>=thisSniffBins(n,1),1,'first');
            x2 = find(TraceTimestamps>=thisSniffBins(n,2),1,'first') - 1;
            OdorMedian(s,n) = median(OdorLocations(x1:x2)); % median location of the odor location in that bin
            OdorStd(s,n) = std(OdorLocations(x1:x2)); % sd of odor location in that bin
        end
    end

    %%
    %figure;
    figure('Name',['unit ',num2str(SelectedUnit)]); %,'Color','none');

    %OdorColors = brewermap(5,'*Paired');
%     MyColors1 = brewermap(8,'YlOrRd');
%     OdorColors(1,:) = Plot_Colors('b'); % ITI sniffs
%     OdorColors(2,:) = MyColors1(4,:); % snifftype = 0 % approach
%     OdorColors(3,:) = MyColors1(6,:); % snifftype = 1 % settle
%     OdorColors(4,:) = MyColors1(7,:); % snifftype = 2 % at target
    
    OdorColors(1,:) = [0.6 0.6 0.6];
    OdorColors(2:4,:) = brewermap(3,'Dark2');


    FRColors = brewermap((1+max(BinnedPSTH(:))),'RdPu');
    %FRColors(1:2,:) = [];

    OdorStates = [-1 1 2 3]; %[-1 0 1 2 3];
    for x = 1:numel(OdorStates)
        WhichSniffs = find(SniffsOEPS(:,6)==OdorStates(x));

        for WhichBin = 1:size(BinnedPSTH,2)
            subplot(numel(OdorStates)+3,size(BinnedPSTH,2),(x-1)*size(BinnedPSTH,2) + WhichBin);
            hold on

            for k = 1:numel(WhichSniffs) % every sniff
                %plot3(SniffsOEPS(WhichSniffs(k),7),SniffsOEPS(WhichSniffs(k),4),BinnedPSTH(WhichSniffs(k),WhichBin),'.','MarkerSize',12,'Color',FRColors(1+BinnedPSTH(WhichSniffs(k),WhichBin),:));
                %plot3(SniffsOEPS(WhichSniffs(k),4),...
                 plot3(OdorMedian(WhichSniffs(k),WhichBin),...
                    SniffsOEPS(WhichSniffs(k),5),...
                    BinnedPSTH(WhichSniffs(k),WhichBin),...
                    '.','MarkerSize',12,'Color',FRColors(1+BinnedPSTH(WhichSniffs(k),WhichBin),:));
            end

            view([0.1 0.2 0.3]);


            % FR vs. sniff duration
            subplot(numel(OdorStates)+3,size(BinnedPSTH,2),(numel(OdorStates))*size(BinnedPSTH,2) + WhichBin);
            hold on

            FRvalues = [];
            for d = 1:8
                thisDurationSniffs = WhichSniffs(find(ceil(SniffsOEPS(WhichSniffs,8)/0.025)==d));
                if ~isempty(thisDurationSniffs)
                    FRvalues(d,1) = mean(BinnedPSTH(thisDurationSniffs,WhichBin));
                    FRvalues(d,2) = std(BinnedPSTH(thisDurationSniffs,WhichBin));
                    FRvalues(d,3) = FRvalues(d,2)/sqrt(numel(thisDurationSniffs));
                    FRvalues(d,4) = mean(SniffsOEPS(thisDurationSniffs,8));
                    FRvalues(d,5) = std(SniffsOEPS(thisDurationSniffs,8));
                else
                    FRvalues(d,1:5) = nan; 
                end
            end
            if ~isempty(FRvalues)
                errorbar(FRvalues(:,4),FRvalues(:,1),FRvalues(:,3),'color',OdorColors(x,:),'LineWidth',1);
            end
            
            % FR vs. sniff slope
            subplot(numel(OdorStates)+3,size(BinnedPSTH,2),(numel(OdorStates)+1)*size(BinnedPSTH,2) + WhichBin);
            hold on

            FRvalues = [];
            for d = 1:8
                thisSlopeSniffs = WhichSniffs(find(ceil(SniffsOEPS(WhichSniffs,5)/2)==d));
                if ~isempty(thisSlopeSniffs)
                    FRvalues(d,1) = mean(BinnedPSTH(thisSlopeSniffs,WhichBin));
                    FRvalues(d,2) = std(BinnedPSTH(thisSlopeSniffs,WhichBin));
                    FRvalues(d,3) = FRvalues(d,2)/sqrt(numel(thisSlopeSniffs));
                    FRvalues(d,4) = mean(SniffsOEPS(thisSlopeSniffs,5));
                    FRvalues(d,5) = std(SniffsOEPS(thisSlopeSniffs,5));
                else
                    FRvalues(d,1:5) = nan; 
                end
            end
            if ~isempty(FRvalues)
                errorbar(FRvalues(:,4),FRvalues(:,1),FRvalues(:,3),'color',OdorColors(x,:),'LineWidth',1);
            end

            % FR vs. odor location
            subplot(numel(OdorStates)+3,size(BinnedPSTH,2),(numel(OdorStates)+2)*size(BinnedPSTH,2) + WhichBin);
            hold on

            FRvalues = [];
            count = 0;
            for d = -10:1:10
                count = count + 1;
                %thisLocationSniffs = WhichSniffs(find(ceil(SniffsOEPS(WhichSniffs,4)/10)==d));
                thisLocationSniffs = WhichSniffs(find(ceil(OdorMedian(WhichSniffs,WhichBin)/10)==d));
                if ~isempty(thisLocationSniffs)
                    FRvalues(count,1) = mean(BinnedPSTH(thisLocationSniffs,WhichBin));
                    FRvalues(count,2) = std(BinnedPSTH(thisLocationSniffs,WhichBin));
                    FRvalues(count,3) = FRvalues(count,2)/sqrt(numel(thisLocationSniffs));

%                     FRvalues(count,4) = mean(SniffsOEPS(thisLocationSniffs,4));
%                     FRvalues(count,5) = std(SniffsOEPS(thisLocationSniffs,4));

                    FRvalues(count,4) = mean(OdorMedian(thisLocationSniffs,WhichBin));
                    FRvalues(count,5) = std(OdorMedian(thisLocationSniffs,WhichBin));
                else
                    FRvalues(count,1:5) = nan; 
                end
            end
            if ~isempty(FRvalues)
                errorbar(FRvalues(:,4),FRvalues(:,1),FRvalues(:,3),'color',OdorColors(x,:),'LineWidth',1);
            end

        end

    end

    set(gcf,'Position',[2399  40  820  1000]);

    if savefigs
        saveas(gcf,fullfile(FigPath,['unit ',num2str(SelectedUnit),'.png']));
        close(gcf);
    end

    disp('done');
end