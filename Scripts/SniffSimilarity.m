% script to take sniffs of particular type
% compute single sniff PSTHs
% compute pairwise distances between the sniff PSTHs

%% load the data
MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q8/Q8_20221204_r0_processed.mat';
MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q9/Q9_20221116_r0_processed.mat';
[TrialAligned, TrialInfo, ...
    ReplayAligned, ReplayInfo, ...
    TuningAligned, TuningInfo, ...
    AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

% get sniff time stamps and info for the sniffs we want to plot
[SelectedSniffs] = SelectSniffs(TrialAligned, TrialInfo, [1 2 3], 'includeITI', 1);

%% sort the sniffs
sortorder = 1; % duration of current sniff
for whichodor = 1:3
    SelectedSniffs{whichodor} = ...
        SortSniffs(SelectedSniffs{whichodor}, sortorder);
end

%% get spike rasters
whichUnit = 9;
for whichodor = 1:3
    % get the spike raster
    [nSniffs{whichodor},AllFR{whichodor},SpikesPlot{whichodor}] = ...
        SniffAlignedPlot(SelectedSniffs{whichodor}, AllUnits.Spikes{whichUnit}, ...
        'plotevents', 0);
end

%% get smoothened psths for individual sniffs
% kernelsize = 10; % in ms
% downsample = 1000; % stick to ms resolution
% taxis = -500:500;  % make a time axis of 1000 ms
% gauss_kernel = normpdf(taxis, 0, kernelsize);
% gauss_kernel = gauss_kernel ./ sum(gauss_kernel);
bufferwindow = 0.1; % seconds
bufferbins = bufferwindow*1000;
myPSTH = [];
binsize = 10; % 20; 50;

for whichodor = 1:3
    % find the max sniff duration for this odor
    maxduration = max(SelectedSniffs{whichodor}(:,15)); % in seconds
    myPSTH{whichodor} = nan(size(SelectedSniffs{whichodor},1),floor(1000*maxduration));
    for n = 1:size(SelectedSniffs{whichodor},1) % every sniff
        sniffduration   = SelectedSniffs{whichodor}(n,15);
        whichspikes     = find(SpikesPlot{whichodor}(:,2) == n);
        thissniffspikes = SpikesPlot{whichodor}(whichspikes,1);
        % truncate spikes
        thissniffspikes(thissniffspikes<-bufferwindow) = [];
        thissniffspikes(thissniffspikes>(sniffduration + bufferwindow)) = [];
        thissniffspikes = thissniffspikes + bufferwindow; % no negative spiketimes
        thissniffpsth   = MakePSTH(thissniffspikes', 0, [], 'kernelsize', binsize);
        
        thissnifflength = floor((sniffduration+bufferwindow)*1000); 
        if length(thissniffpsth) < thissnifflength
            thissniffpsth = [thissniffpsth, zeros(1,thissnifflength-length(thissniffpsth))];
        end
        myPSTH{whichodor}(n,1:(thissnifflength-bufferbins)) = thissniffpsth((1+bufferbins):thissnifflength);
        PSTHLength{whichodor}(n,1) = thissnifflength-bufferbins;
    end
end

%% compute similarity for only the ITI sniffs
for whichodor = 1:3
    whichsniffs = 1:size(SelectedSniffs{whichodor},1); %find(SelectedSniffs{whichodor}(:,4)==-1);
    Pairs = nchoosek(whichsniffs,2);
    for p = 1:size(Pairs,1)
        snifflength = min(PSTHLength{whichodor}(Pairs(p,1:2),1));
        PSTHs = myPSTH{whichodor}(Pairs(p,1:2),1:snifflength)';
        c = corrcoef(PSTHs);
        Pairs(p,3) = c(1,2); % correlation
        Pairs(p,4) = sqrt(mean(diff(PSTHs').^2)); % rmse
        
        Pairs(p,5) = diff(SelectedSniffs{whichodor}(Pairs(p,1:2),15)); % difference in duration
        Pairs(p,6) = diff(SelectedSniffs{whichodor}(Pairs(p,1:2),17)); % difference in previous sniff duration
        Pairs(p,7) = sum(SelectedSniffs{whichodor}(Pairs(p,1:2),2)); % difference in odor state
        Pairs(p,8) = diff(SelectedSniffs{whichodor}(Pairs(p,1:2),13)); % difference in odor location
        Pairs(p,9) = diff(SelectedSniffs{whichodor}(Pairs(p,1:2),12)); % difference in previous sniff odor location
    end
    Distances{whichodor} = Pairs;
end

%% plotting

for fignum = 1:3
    figure;
    colormap(brewermap([100],'*Blues'));
    switch fignum
        case 1
            mycol = 8; % difference in odor location
            set(gcf, 'Name', 'difference in odor location');
            XYBins = [25 25];
        case 2
            mycol = 5; % difference in sniff duration
            set(gcf, 'Name', 'difference in sniff duration');
            XYBins = [25 25];
        case 3
            mycol = 6; % difference in previous sniff duration
            set(gcf, 'Name', 'difference in previous sniff duration');
            XYBins = [25 25];
    end
    
    for whichodor = 1:3
        for pairtype = 1:3
            switch pairtype
                case 1
                    % ITI sniffs
                    f = find(Distances{whichodor}(:,7)==-2);
                    subplot(3,3,whichodor);
                case 2
                    % in Trial sniffs
                    f = find(Distances{whichodor}(:,7)==2*whichodor);
                    subplot(3,3,whichodor+3);
                case 3
                    % mixed pairs
                    f = find(Distances{whichodor}(:,7)==(whichodor-1));
                    subplot(3,3,whichodor+6);
            end
            N = histcounts2(abs(Distances{whichodor}(f,mycol)),Distances{whichodor}(f,4),XYBins);
            N = N./sum(N,2);
            %imagesc(flipud(N'));
            imagesc(flipud(N'),[0 .15]);
        end
    end
end
