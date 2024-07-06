%% script to extract the sniff aligned closed-loop spiking responses
%  of a given neuron and fit the ITI, Air and Odor kernels

%% Step 1: Get the spiking data and sniff parameters

%SessionName = 'O9_20220630_r0'; MyUnits = [];
%SessionName = 'O8_20220702_r0'; MyUnits = [];
%SessionName = 'S12_20230804_r0'; MyUnits = [];
%SessionName = 'Q4_20221109_r0'; MyUnits = [];
%SessionName = 'Q8_20221204_r0'; MyUnits = [];
%SessionName = 'Q9_20221116_r0'; MyUnits = [];
%SessionName = 'Q3_20221019_r0'; MyUnits = [];
%SessionName = 'Q4_20221109_r0'; MyUnits = [2 55 12 69 4 9 19 14 16 26 41 10]; %Q4
SessionName = 'S12_20230731_r0'; MyUnits = [2 3 6 9 10 11 16 17 18 19 25 26 28 29 31 34 39 42 44 46 47 48 49 50 54 58 69 70 72 73 74 85 87 88 95 97]; % S12
%SessionName = 'Q9_20221116_r0'; MyUnits = [1 11 15 18 19 23 28 29 36 39 43 49 94]; %Q9
%SessionName = 'Q8_20221204_r0'; MyUnits = [49 54 104]; %Q8
%SessionName = 'Q3_20221019_r0'; MyUnits = [16 18]; %Q3
%SessionName = 'O3_20210927_r0'; MyUnits = [2 3 7 9 11 13 14 15 19 25 27 28 33 35 42 43 45 48 54 58 59 62 66 72] % O3

MySession = fullfile('/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/', ...
                        SessionName(1:regexp(SessionName,'_','once')-1), ...
                        [SessionName,'_processed.mat']);
FigPath = ['/home/priyanka/Desktop/sniffPSTHPredictions/', SessionName(1:regexp(SessionName,'_','once')-1)];
if ~exist(FigPath,'dir')
    mkdir(FigPath);
    fileattrib(FigPath, '+w','a');
end

[TrialAligned, TrialInfo, ...
    ReplayAligned, ReplayInfo, ...
    TuningAligned, TuningInfo, ...
    AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

%% make continous time vectors
AllSniffs = [];
for trial = 1:size(TrialAligned.Sniffs,2)
    thisTrialSniffs = TrialAligned.Sniffs{trial}(:,[7 8 2 13]); % inhalationStart inhalationEnd Odor/AirState MotorLocation
    % add the actual trial start time
    thisTrialSniffs(:,1:2) = thisTrialSniffs(:,1:2) + TrialAligned.RefTS(trial);
    AllSniffs = vertcat(AllSniffs,thisTrialSniffs);
end
% everything above is in OEPS timebase

% remove duplicates if any
if size(unique(AllSniffs(:,1:2),'rows'),1) < size(AllSniffs,1)
    keyboard;
    U = unique(AllSniffs(:,1:2),'rows');
    AllSniffs = AllSniffs(U,:);
end

% make a long vector @10ms resolution until the last sniff end + 1s
PSTHbinsize = 10; 
deltasniffs = 1; % make all sniffs the same duration

nBins = ceil(AllSniffs(end,2)*1000/PSTHbinsize);
InputVector = zeros(6,nBins); % sniff, air, odor1, odor2, odor3, motor location

for s = 1:size(AllSniffs,1) % every sniff
    ts = AllSniffs(s,1:2); % inhalation start and end in seconds 
    ts = ts*1000/PSTHbinsize; % in bins
    ts = round(ts);
    if deltasniffs
        ts(2) = ts(1);
        %ts(2) = ts(1) + (50/PSTHbinsize) - 1;
    end
    InputVector(1,ts(1):ts(2)) = 1; % inhalation period
    if AllSniffs(s,3) >= 0
        InputVector(2,ts(1):ts(2)) = 1; % Sniffs with Manifold Air ON 
        InputVector((2 + AllSniffs(s,3)),ts(1):ts(2)) = 1; % only pull up inhalations in the corresponding odor vector
    end    
    InputVector(6,ts(1):ts(2)) = AllSniffs(s,4); % keep track of odor location in each sniff
end

%% Get the PSTHs at the same resolution
PSTHs = zeros(size(AllUnits.Spikes,2),nBins);
for n = 1:size(AllUnits.Spikes,2) % every unit
    thisUnitSpikeTimes = AllUnits.Spikes{n}; % in seconds
    thisUnitSpikeTimes = floor(thisUnitSpikeTimes*1000/PSTHbinsize) + 1; % in binIDs
    % ignore spiketimes larger than the largest sniff time
    thisUnitSpikeTimes(thisUnitSpikeTimes>nBins,:) = [];
    [C,~,ic] = unique(thisUnitSpikeTimes);
    bin_counts = accumarray(ic,1);
    if ~isempty(C)
        PSTHs(n,C) = bin_counts;
        %PSTHs(n,C) = 1000*bin_counts/bindownto; % FR in Hz
    end
end

%% set up the model
% start with random kernels
unitid = 58;
kernelLength = 700; % in ms
StartingKernels = zeros(1,(5*(kernelLength/PSTHbinsize)) + 2); % 5 kernels + baseline + locationcoeff
[kernelsout_cont{unitid,1}] = GetSniffKernelsContinuous(StartingKernels, InputVector, PSTHs(unitid,:), 'binsize', PSTHbinsize, 'rectifyFR', 0);

%% get the predicted PSTH
[baseline,kernels,locationcoef] = ParseSniffKernels(kernelsout_cont{unitid,1});
[zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,InputVector);

%% get the ITI PSTH
In_On = find(diff(InputVector(1,:))>0);
In_On(2,:) = [In_On(1,2:end) NaN];
In_On(:,find(InputVector(2,In_On(1,:)+1)>0)) = [];
ITIPSTH = [];
maxbins = max(diff(In_On,1)) + 1 + 20;
for x = 1:size(In_On,2)
    bin1 = In_On(1,x) - 9;
    bin2 = bin1 + maxbins;
    myPSTH = zdata(1,bin1:bin2);
    bins2nan = In_On(2,x) - In_On(1,x) + 20;
    myPSTH(1,bins2nan:end) = NaN;
    ITIPSTH = vertcat(ITIPSTH, myPSTH);
end

[~,I] = sort(diff(In_On,1));
ITIPSTH = ITIPSTH(I,:);

tokeep = sort(randperm(4100,2000));
ITIPSTH = ITIPSTH(tokeep,:);

%%
dt = 0.010; % in seconds
SpikesPlot = [];
for x = 1:size(ITIPSTH,1)   
    % get the predicted spike times
    thisSniffSpikes = PechePourPoisson(100*ITIPSTH(x,:),dt) - 0.1;
    SpikesPlot = vertcat(SpikesPlot, [thisSniffSpikes x*ones(numel(thisSniffSpikes),1)]);
    %plot(T, x + (0*T), 'k', 'MarkerSize', 2);
end


% plot spikes
plot(SpikesPlot(:,1),SpikesPlot(:,2),'.k','Markersize', 2);
%set(gca,'YLim', [0 2200], 'XLim',[-0.05 0.75],'TickDir','out');
set(gca,'XLim',[-0.05 0.75],'TickDir','out');

%%
% get sniff time stamps and info for the sniffs we want to plot
[SelectedSniffs] = SelectSniffs_forKernelFits(TrialAligned, TrialInfo, [1 2 3]);

if isempty(MyUnits)
    MyUnits = 1:size(AllUnits.ChannelInfo,1); % all Units 
end

%%
PSTHbinsize = 10;
for unitcount = 1:numel(MyUnits)
    whichUnit = MyUnits(unitcount);
    
    SniffPSTHs = []; SniffParams = [];
    
    % get spike rasters
    for whichodor = 1:3
        % get the spike raster
        [SniffPSTHs_temp,SniffParams_temp] = ...
            SniffAlignedFRs(SelectedSniffs{whichodor}, AllUnits.Spikes{whichUnit}, ...
            'psthtype', 'binned', 'binsize', PSTHbinsize);
        
        SniffDim = size(SniffPSTHs_temp);
        SniffPSTHs(end+(1:SniffDim(1)),1:SniffDim(2)) = SniffPSTHs_temp;
        SniffParams = vertcat(SniffParams, SniffParams_temp);
    end
    % sniff params
    % 1-4: [currsniffstate currsniffloc currinhend currsniffend ...
    % 5-6:  currsniffTrialID currsniffIndex ...
    % 7:10: prevsniffstate prevsniffloc previnhstart previnhend]
    
    % Split sniffs into 2 halves
    nsniffs = size(SniffParams,1);
    AllsniffIDs = randperm(nsniffs);
    sniffIDs{1} = AllsniffIDs(1:floor(nsniffs/2));
    sniffIDs{2} = AllsniffIDs((1+floor(nsniffs/2)):end);
    
    %SniffsUsed{unitcount} = sniffIDs;
    Sniffs.Params{unitcount} = SniffParams;
    Sniffs.PSTHs{unitcount}  = SniffPSTHs;
    
    for sniffset = 1:2
      
        Sniffs.IDs{unitcount}{sniffset} = sniffIDs{sniffset};
              
        SniffPSTHs_  = SniffPSTHs(sniffIDs{sniffset},:);
        SniffPSTHs_  = SniffPSTHs_(:,1:(1+max(SniffPSTHs_(:,1))));
        SniffParams_ = SniffParams(sniffIDs{sniffset},:);
        
        % starting kernels
        [StartingKernels] = InitialKernelEstimates(SniffPSTHs_, SniffParams_, ...
            'kernellength', 700, 'binsize', PSTHbinsize);
        
        %kernelsIn{unitcount} = StartingKernels;
        [kernelsout{unitcount,sniffset},resnorm{unitcount,sniffset},residual{unitcount,sniffset},exitflag,output] = ...
            GetSniffKernels(StartingKernels, SniffParams_, SniffPSTHs_, ...
            'binsize', PSTHbinsize, 'rectifyFR', 0);
        
    end
end

%%
SavePath = fullfile(FigPath,[SessionName(1:regexp(SessionName,'_r','once')-1),'_sniffs.mat']);
save(SavePath,'SessionName','MyUnits','Sniffs','kernelsout','PSTHbinsize');

%% compare PSTH
%MyUnits = sort(MyUnits);
for unitcount = 1:numel(MyUnits)
    whichUnit = MyUnits(unitcount);
    
    SniffPSTHs = []; SniffParams = [];
    
    % get spike rasters
    for whichodor = 1:3
        % get the spike raster
        [SniffPSTHs_temp,SniffParams_temp] = ...
            SniffAlignedFRs(SelectedSniffs{whichodor}, AllUnits.Spikes{whichUnit}, ...
            'psthtype', 'binned', 'binsize', PSTHbinsize);
        
        SniffDim = size(SniffPSTHs_temp);
        SniffPSTHs(end+(1:SniffDim(1)),1:SniffDim(2)) = SniffPSTHs_temp;
        SniffParams = vertcat(SniffParams, SniffParams_temp);
    end
    
    % predicted PSTHs
    [PSTHOut] = GetPredictedSniffPSTH(SniffParams,kernelsout{unitcount},'binsize', PSTHbinsize);
    
    %PSTHOut = PSTHOut*10;
    
    figure('Name',['unit ',num2str(whichUnit)],'Color','none');
    colormap(brewermap([100],'*YlGnBu'));
    odorstates = [-1 1 2 3]; 
    
    % plot stuff - sort sniffs by duration
    alims = [0 0]; blims = [0 0];
    for i = 1:4
        whichsniffs = find(SniffParams(:,1)==odorstates(i));
        
        % sort sniffs by duration
        [~,s] = sort(SniffParams(whichsniffs,4));
        
        subplot(2,8,(i*2)-1);
        imagesc(flipud(SniffPSTHs(whichsniffs(s),2:end)));
        alims = max([alims; get(gca,'CLim')]);
        
        subplot(2,8,(i*2)-0);
        imagesc(flipud(PSTHOut(whichsniffs(s),1:end)));
        blims = max([blims; get(gca,'CLim')]);
        
        % sort sniffs by odor location
        [~,s] = sort(SniffParams(whichsniffs,2));
        
        subplot(2,8,8 +(i*2)-1);
        imagesc(flipud(SniffPSTHs(whichsniffs(s),2:end)));
        
        subplot(2,8,8 + (i*2)-0);
        imagesc(flipud(PSTHOut(whichsniffs(s),1:end)));
    end
    for i = 2:2:16
        subplot(2,8,i-1);
        set(gca,'Clim',[0 alims(2)]);
        subplot(2,8,i);
        set(gca,'Clim',[0 blims(2)]);
    end
    set(gcf,'Position',[143 403 1653 420]);
    saveas(gcf,fullfile(FigPath,['unit ',num2str(whichUnit),'.png']));
    close all
end

%%
cf = []; 
for i = 1:numel(MyUnits) 
    cf(i,:) = kernelsout{i}([1 end]); 
    %cf(i,:) = kernelsout{i}([1 end-2:end]); 
end
figure('Name','location-coeffs & baselines')
subplot(2,1,1);
%plot(cf(:,2:4),'-o'); %,'XTicklabel',MyUnits);
plot(cf(:,2),'-o'); 
xticks(1:numel(MyUnits));
xticklabels(num2str(MyUnits'));

subplot(2,1,2);
plot(cf(:,1),'-o'); %,'XTicklabel',MyUnits);
xticks(1:numel(MyUnits));
xticklabels(num2str(MyUnits'));
saveas(gcf,fullfile(FigPath,['fittedcoeffs.png']));

close all

disp('done!');

