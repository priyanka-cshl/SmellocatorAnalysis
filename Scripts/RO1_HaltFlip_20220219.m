%% for O3
SessionPath = 'O3/O3_20210927_r0_processed.mat';
ChosenUnits = [15 17 18]; % 28 63 56]; %MyUnits = [8 35 28 55 39]; 
ChosenUnits = [15 17 18 14 25 54];

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

%% Get Odor response values from TrialStart and from pertubation period

TrialInfo.TargetEntry = NaN*TrialInfo.Odor;
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);


whichOdor = 1;
AlignTo = 6; % perturbation start
MyColors1 = brewermap(15,'Greys');
H = figure;
nCols = numel(ChosenUnits)+2;
for i = 1:numel(ChosenUnits)
    whichUnit = ChosenUnits(i);
    subplot(3,nCols,i); hold on
    [myFRs, BinOffset, whichTZ] = PlotHaltFlips(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo);
    set(gca, 'XLim', [-1.2 6],'TickDir','out','XTick',[]);
    title(['unit# ',num2str(whichUnit)]);
    
    subplot(3,nCols,i+nCols); hold on
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
    
    subplot(3,nCols,i+2*nCols); hold on
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
    % get FRs aligned to odor start as well - 
    [~, FRs, BinOffset] = PlotFullSession(-whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, 2);
    for t = 1:size(FRs,1)
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(t,:),'Color',MyColors1(t,:),'Linewidth',1);
    end
    set(gca, 'XLim', [-1.2 6],'TickDir','out');
    
    % add odor start response to previous plot as well
    subplot(3,nCols,i+nCols); hold on
    plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(whichTZ,:),'color',Plot_Colors('r'),'Linewidth',1);
    set(gca, 'XLim', [-1.2 6],'TickDir','out','XTick',[]);
end

%% behavior plot
perturbationTrials = find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip'));
perturbationTrials = [perturbationTrials TrialInfo.TargetZoneType(perturbationTrials) Events(perturbationTrials,5)];
perturbationTrials = sortrows(perturbationTrials,2);
nZones = unique(perturbationTrials(:,2));

H1 = figure;
for i = 1:numel(nZones)
    figure(H1);
    whichTZ = nZones(i);
    subplot(4,numel(nZones),i); hold on
    % plot all perturbation trials, aligned to perturbation start
    thisZoneTrace = [];
    thisZoneTrials = find(perturbationTrials(:,2)==whichTZ);
    for j = 1:numel(thisZoneTrials)
        foo = Traces.Lever{perturbationTrials(thisZoneTrials(j),1)};
        % delete samples such as to align to perturbation start
        nSamps = round(perturbationTrials(thisZoneTrials(j),3)*SampleRate);
        foo(1:nSamps,:) = [];
        plot(foo,'color',Plot_Colors('t'));
        thisZoneTrace(1:length(foo),j) = foo; 
    end
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out','XTick',[]);
    subplot(4,numel(nZones),i+numel(nZones)); hold on
    MyShadedErrorBar(1:size(thisZoneTrace,1),mean(thisZoneTrace,2)',std(thisZoneTrace'),Plot_Colors('t'),{},0.5);
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out','XTick',[]);
    
    % control trials
    controlTrials = intersect(find(TrialInfo.TargetZoneType==whichTZ),...
                              find(TrialInfo.Odor == whichOdor));
    controlTrials = intersect(find(~strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')), ...
                              controlTrials);
    controlTrials(:,3) = Events(controlTrials,5);
    ControlTrace = [];
    subplot(4,numel(nZones),i+2*numel(nZones)); hold on
    for j = 1:size(controlTrials,1)
        foo = Traces.Lever{controlTrials(j,1)};
        % delete samples such as to align to perturbation start
        nSamps = round(controlTrials(j,3)*SampleRate);
        foo(1:nSamps,:) = [];
        plot(foo,'color',Plot_Colors('k'));
        ControlTrace(1:length(foo),j) = foo; 
    end
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out','XTick',[]);
    subplot(4,numel(nZones),i+3*numel(nZones)); hold on
    MyShadedErrorBar(1:size(ControlTrace,1),mean(ControlTrace,2)',std(ControlTrace'),Plot_Colors('k'),{},0.5);
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out');
    
    if i == 3
        % for this zone lets get concatenated trials
        for j = 1:numel(controlTrials(:,1))
            Traces.HaltFlip(controlTrials(j,1)) = {0*Traces.Trial{controlTrials(j,1)}};
            Traces.TargetZone(controlTrials(j,1)) = {TrialInfo.TargetZoneType(controlTrials(j,1)) + ...
                0*Traces.Trial{controlTrials(j,1)}};
        end
        for j = 1:numel(thisZoneTrials(:,1))
            foo = 0*Traces.Trial{perturbationTrials(thisZoneTrials(j),1)};
            idx = TrialInfo.Perturbation{perturbationTrials(thisZoneTrials(j),1),2};
            idx(:,1:2) = idx(:,1:2) + SampleRate*startoffset;
            foo(idx(1):idx(2)) = -idx(3);
            Traces.HaltFlip(perturbationTrials(thisZoneTrials(j),1)) = {foo};
            Traces.TargetZone(perturbationTrials(thisZoneTrials(j),1)) = ...
                {TrialInfo.TargetZoneType(perturbationTrials(thisZoneTrials(j),1)) + ...
                0*foo};
        end
        allTrials = [controlTrials(:,1); perturbationTrials(thisZoneTrials,1)];
        
        B1 = figure;
        figure(B1);
        [TracesOut] = ConcatenateTraces(Traces, allTrials, SampleRate*startoffset);
        timestamps = (1:size(TracesOut.Lever{1},1))'/SampleRate;
        TZ = [TracesOut.TargetZone{1} TracesOut.TargetZone{1}];
        TZ = [TargetZones(TZ(:,1),1) TargetZones(TZ(:,1),3)];
        Trial = TracesOut.Trial{1};
        Trial(Trial~=whichOdor) = 0;
        PlotBehaviorPerturbations(timestamps,TracesOut.Lever{1},...
            TracesOut.Sniffs{1},TracesOut.Licks{1},TracesOut.Rewards{1},...
            Trial,...
            TZ, ...
            TracesOut.HaltFlip{1},5);
        
    end
end


%%
for tz = 1:12
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
                    find(TrialInfo.Odor==1));
    whichTrials = intersect(whichTrials,...
                    find(TrialInfo.TargetZoneType==tz));
    OdorONPeriod(tz,1) = median(TrialInfo.OdorStart(whichTrials,1));
    foo = cellfun(@(x) median(x(491:500)), Traces.Motor(whichTrials),'UniformOutput', false);
    LocationTrialStart(tz,1) = round(median(cell2mat(foo)));
end
% convert to Bins
OdorONPeriod = ceil(OdorONPeriod*SampleRate); 

% Align to:
% 1 - Trial ON, 2 - Odor ON, 3 - Trial OFF, 4 - Reward, 5 - 1st TZ entry, 
% 6 - Perturbation Start
%%
mywin = 1000; % in ms - to be defined by user
stepsize = 1000/SampleRate; % method used to account for potential sampledrop
for i = 1:size(SingleUnits,2)
    % get spike counts aligned to odor ON
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 2);
    % count spikes in a 1000 ms (mywin) bin after odor ON
    bins = -BinOffset + [stepsize:stepsize:mywin]; %in ms
    bins = bins/stepsize; % to convert to indices
    for tz = 1:12
        % get odor response during 'windowsize' after odor ON
        AreaUnderCurve(tz,1,i) = sum(AlignedFRs(tz,bins));
    end
    % get spike counts aligned to perturbation start
    [~, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = ... 
        PlotFullSession(-i, whichOdor, AlignedSpikes, Events, TrialInfo, 6);
    % count spikes in a 1000 ms bin after odor ON
    bins = -BinOffset + [stepsize:stepsize:mywin]; %in ms
    bins = bins/stepsize; % to convert to indices
    for tz = 1:12
        AreaUnderCurve(tz,2,i) = sum(AlignedPerturbationFRs(tz,bins));
    end
end
AreaUnderCurve = AreaUnderCurve/(mywin/stepsize);

%% Plot the effect
figure;
subplot(1,2,1); 
hold on
axis square
subplot(1,2,2); 
hold on
axis square
for i = 1:size(SingleUnits,2)
    % mean perturbation response vs. mean mirror location response 
    % = roughly LocationTrialStart(11:12)
    subplot(1,2,1);
    plot(mean(AreaUnderCurve(11:12,1,i)),mean(AreaUnderCurve([1 5 9],2,i)),...
        'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
    % halt perturbations happen for tz type 1, 5 and 9
    
    line([0 40],[0 40],'Color','k')
    subplot(1,2,2);
    %comparing to halt response to max response at this location no matter
    %what is the trial type
    plot(max(AreaUnderCurve(:,1,i)),mean(AreaUnderCurve([1 5 9],2,i)),...
        'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
    line([0 40],[0 40],'Color','k')
end


%% add PCX4 examples to the plot
SessionPath = 'PCX4/PCX4_20210713_r0_processed.mat';
PCXUnits = [7 16]; % 28 63 56]; %MyUnits = [8 35 28 55 39]; 

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

%% Get Odor response values from TrialStart and from pertubation period

TrialInfo.TargetEntry = NaN*TrialInfo.Odor;
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

figure(H);
for i = numel(ChosenUnits)+(1:numel(PCXUnits))
    whichUnit = PCXUnits(i - numel(ChosenUnits));
    subplot(3,nCols,i); hold on
    [myFRs, BinOffset, whichTZ] = PlotHaltFlips(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo);
    set(gca, 'XLim', [-1.2 6],'TickDir','out','XTick',[]);
    title(['unit# ',num2str(whichUnit)]);
    
    subplot(3,nCols,i+nCols); hold on
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
    
    subplot(3,nCols,i+2*nCols); hold on
    plot((1:size(myFRs,2))*0.002+BinOffset/1000,myFRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
    % get FRs aligned to odor start as well - 
    [~, FRs, BinOffset] = PlotFullSession(-whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, 2);
    for t = 1:size(FRs,1)
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(t,:),'Color',MyColors1(t,:),'Linewidth',1);
    end
    set(gca, 'XLim', [-1.2 6],'TickDir','out');
    
    % add odor start response to previous plot as well
    subplot(3,nCols,i+nCols); hold on
    plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(whichTZ,:),'color',Plot_Colors('r'),'Linewidth',1);
    set(gca, 'XLim', [-1.2 6],'TickDir','out','XTick',[]);
end


%% behavior plot
perturbationTrials = find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip'));
perturbationTrials = [perturbationTrials TrialInfo.TargetZoneType(perturbationTrials) Events(perturbationTrials,5)];
perturbationTrials = sortrows(perturbationTrials,2);
nZones = unique(perturbationTrials(:,2));

H2 = figure
for i = 1:numel(nZones)
    figure(H2);
    whichTZ = nZones(i);
    subplot(4,numel(nZones),i); hold on
    % plot all perturbation trials, aligned to perturbation start
    thisZoneTrace = [];
    thisZoneTrials = find(perturbationTrials(:,2)==whichTZ);
    for j = 1:numel(thisZoneTrials)
        foo = Traces.Lever{perturbationTrials(thisZoneTrials(j),1)};
        % delete samples such as to align to perturbation start
        nSamps = round(perturbationTrials(thisZoneTrials(j),3)*SampleRate);
        foo(1:nSamps,:) = [];
        plot(foo,'color',Plot_Colors('t'));
        thisZoneTrace(1:length(foo),j) = foo; 
    end
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out','XTick',[]);
    subplot(4,numel(nZones),i+numel(nZones)); hold on
    MyShadedErrorBar(1:size(thisZoneTrace,1),mean(thisZoneTrace,2)',std(thisZoneTrace'),Plot_Colors('t'),{},0.5);
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out','XTick',[]);
    
    % control trials
    controlTrials = intersect(find(TrialInfo.TargetZoneType==whichTZ),...
                              find(TrialInfo.Odor == whichOdor));
    controlTrials = intersect(find(~strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')), ...
                              controlTrials);
    controlTrials(:,3) = Events(controlTrials,5);
    ControlTrace = [];
    subplot(4,numel(nZones),i+2*numel(nZones)); hold on
    for j = 1:size(controlTrials,1)
        foo = Traces.Lever{controlTrials(j,1)};
        % delete samples such as to align to perturbation start
        nSamps = round(controlTrials(j,3)*SampleRate);
        foo(1:nSamps,:) = [];
        plot(foo,'color',Plot_Colors('k'));
        ControlTrace(1:length(foo),j) = foo; 
    end
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out','XTick',[]);
    subplot(4,numel(nZones),i+3*numel(nZones)); hold on
    MyShadedErrorBar(1:size(ControlTrace,1),mean(ControlTrace,2)',std(ControlTrace'),Plot_Colors('k'),{},0.5);
    set(gca,'XLim',[0 2000],'YLim',[0 5.5],'TickDir','out');
    
    if i == 3
        % for this zone lets get concatenated trials
        for j = 1:numel(controlTrials(:,1))
            Traces.HaltFlip(controlTrials(j,1)) = {0*Traces.Trial{controlTrials(j,1)}};
            Traces.TargetZone(controlTrials(j,1)) = {TrialInfo.TargetZoneType(controlTrials(j,1)) + ...
                0*Traces.Trial{controlTrials(j,1)}};
        end
        for j = 1:numel(thisZoneTrials(:,1))
            foo = 0*Traces.Trial{perturbationTrials(thisZoneTrials(j),1)};
            idx = TrialInfo.Perturbation{perturbationTrials(thisZoneTrials(j),1),2};
            idx(:,1:2) = idx(:,1:2) + SampleRate*startoffset;
            foo(idx(1):idx(2)) = idx(3);
            Traces.HaltFlip(perturbationTrials(thisZoneTrials(j),1)) = {foo};
            Traces.TargetZone(perturbationTrials(thisZoneTrials(j),1)) = ...
                {TrialInfo.TargetZoneType(perturbationTrials(thisZoneTrials(j),1)) + ...
                0*foo};
        end
        allTrials = [controlTrials(:,1); perturbationTrials(thisZoneTrials,1)];
        
        B2 = figure;
        figure(B2);
        [TracesOut] = ConcatenateTraces(Traces, allTrials, SampleRate*startoffset);
        
        timestamps = (1:size(TracesOut.Lever{1},1))'/SampleRate;
        TZ = [TracesOut.TargetZone{1} TracesOut.TargetZone{1}];
        TZ = [TargetZones(TZ(:,1),1) TargetZones(TZ(:,1),3)];
        Trial = TracesOut.Trial{1};
        Trial(Trial~=whichOdor) = 0;
        PlotBehaviorPerturbations(timestamps,TracesOut.Lever{1},...
            TracesOut.Sniffs{1},TracesOut.Licks{1},TracesOut.Rewards{1},...
            Trial,...
            TZ, ...
            TracesOut.HaltFlip{1},5);
        
    end
end


