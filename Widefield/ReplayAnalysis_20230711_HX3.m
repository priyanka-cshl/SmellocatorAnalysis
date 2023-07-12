%% paths
datapath = '/mnt/data/Behavior/';
imagingpath = '/mnt/data/Widefield/HX3/20230505_r0';
MySession = fullfile(datapath,'HX3','HX3_20230505_r0_processed.mat');

%% get the processed data loaded and assemble open loop traces
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags','WF_timestamps');
           
OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

%% get the corresponding frame ids to use using behavior timestamps

% template
template_TS = OpenLoop.TemplateTraces.Timestamps{1};
% find corresponding frame indices
[~,f1] = min(abs(WF_timestamps.Behavior(:,1)-template_TS(1)));
[~,f2] = min(abs(WF_timestamps.Behavior(:,1)-template_TS(end)));
frame_TS = WF_timestamps.Behavior(f1:f2,1);
frame_idx = f1:f2;
valid_frames = find(frame_idx>0);
nan_locs(:,1) = find(diff(isnan(template_TS))==1) + 1;
nan_locs(:,2) = find(diff(isnan(template_TS))==-1);
for i = 1:size(nan_locs,1)
    idx1 = find(frame_TS > template_TS(nan_locs(i,1)-1),1,'first');
    idx2 = find(frame_TS < template_TS(nan_locs(i,2)+1),1,'last');
    frame_idx(idx1:idx2) = -frame_idx(idx1:idx2);
end

% active replays
for r = 1:size(OpenLoop.ReplayTraces.Timestamps{1},2)
    TS(r,1) = find(OpenLoop.ReplayTraces.Timestamps{1}(:,r)~=0,1,'first'); % traces are aligned to replay end, so beginning can be zeros
    TS(r,2) = find(OpenLoop.ReplayTraces.Timestamps{1}(:,r)~=0,1,'last');
    TS(r,:) = OpenLoop.ReplayTraces.Timestamps{1}(TS(r,:),r);
    
    % find corresponding frame indices
    [~,f1] = min(abs(WF_timestamps.Behavior(:,1)-TS(r,1)));
    [~,f2] = min(abs(WF_timestamps.Behavior(:,1)-TS(r,2)));
    Frames_idx(r,:) = [f1 f2];
end

% passive replays
for p = 1:size(PassiveReplayTraces.Timestamps,2)
    TS(r+p,:) = PassiveReplayTraces.Timestamps{p}([1 end]);
    
    % find corresponding frame indices
    [~,f1] = min(abs(WF_timestamps.Passive(:,1)-TS(r+p,1)));
    [~,f2] = min(abs(WF_timestamps.Passive(:,1)-TS(r+p,2)));
    Frames_idx(r+p,:) = [f1 f2];
end

%% load imaging data
% only selected ROIs
load(fullfile(imagingpath,'selectedROIs.mat'), 'WhichROIs','refPix','C','stackdims');   
% also the whole stack
[WF_stack] = ViewProcessedBinaryStack(imagingpath);

%% plotting ROIs
figure;
imagesc(mean(WF_stack,3));
colormap(brewermap([50],'RdBu'));
hold on
% check locations
for i = 1:length(WhichROIs)
    %plot(WhichROIs(i,2),WhichROIs(i,1),'.k');
    text(WhichROIs(i,2),WhichROIs(i,1),num2str(i));
end

%% extracting replay responses of any given ROI
figure; 
Lims = [30 40];
ROIidx = 8; %8; %5; % 9

% Lims = [40 50];
% ROIidx = 4; %4; %7; % 9
% % % 
Lims = [20 30];
ROIidx = 14; %7; % 9

% template
subplot(3,1,1);
ProcessOpenLoopBasic(OpenLoop, SampleRate, TargetZones, 'plottrials', 1, 'TrialHeight',Lims);
hold on
valid_frames = find(frame_idx>0);
ts = frame_TS;
plot(ts(valid_frames)-ts(1),refPix(frame_idx(valid_frames),ROIidx));

% plot all active replay traces
subplot(3,1,2);
ProcessOpenLoopBasic(OpenLoop, SampleRate, TargetZones, 'plottrials', 1, 'TrialHeight',Lims);
hold on
for i = 1:r
    indices = Frames_idx(i,:);
    ts = WF_timestamps.Behavior(indices(1):indices(2),1);
    plot(ts-ts(1),refPix(indices(1):indices(2),ROIidx));
end

% passive replays
subplot(3,1,3); 
ProcessOpenLoopBasic(OpenLoop, SampleRate, TargetZones, 'plottrials', 1, 'TrialHeight',Lims);
hold on
for i = 1:p
    indices = Frames_idx(i+r,:);
    ts = WF_timestamps.Passive(indices(1):indices(2),1);
    indices = WF_timestamps.Passive(Frames_idx(i+r,:),2);
    plot(ts-ts(1),refPix(indices(1):indices(2),ROIidx));
end

%% convert WhichROIs into linear pix identity
for x = 1:length(WhichROIs)
    ROI_idx(x) = stackdims(1)*(WhichROIs(x,2)-1) + WhichROIs(x,1);
end

figure; 
%% all pixel correlation for replays
nPixels = stackdims(1)*stackdims(2);
nFrames = stackdims(3);
WF_stack = reshape(WF_stack,nPixels,nFrames); % each col is a frame [cols are stacked on top of each other]
Cmat = zeros(nPixels,r+p);
for pix = 1:nPixels
    clear MyTraces nSamps
    CL_trace = WF_stack(pix,frame_idx(valid_frames));
    nSamps(1) = length(CL_trace);
    MyTraces(:,1) = CL_trace;
for i = 1:r
    indices = Frames_idx(i,:);
    Replay_trace = WF_stack(pix,indices(1):indices(2));
    nSamps(i+1) = length(Replay_trace);
    MyTraces(1:length(Replay_trace),i+1) = Replay_trace;
%     nSamps = min(length(CL_trace),length(Replay_trace));
%     
%     Cmat(pix,i) = corr(CL_trace(1:nSamps)',Replay_trace(1:nSamps)');
end
for i = 1:p
    indices = WF_timestamps.Passive(Frames_idx(i+r,:),2);
    Replay_trace = WF_stack(pix,indices(1):indices(2));
    nSamps(i+r+1) = length(Replay_trace);
    MyTraces(1:length(Replay_trace),i+r+1) = Replay_trace;
%     nSamps = min(length(CL_trace),length(Replay_trace));
%     Cmat(pix,i+r) = corr(CL_trace(1:nSamps)',Replay_trace(1:nSamps)');
end

MyTraces(1+min(nSamps):end,:) = [];

% now calculate correlations
MyCorr = corrcoef(MyTraces);
m = size(MyCorr,2) - 1; % no. of comparisons, excluding self    
% correlation of replay to closed loop
Cmat(pix,1:m) = MyCorr(1,2:end);

% between replays
MyCorr(1,:) = [];
MyCorr(:,1) = [];

OL_corrs = MyCorr(1:r,1:r); % only the active replays
PR_corrs = MyCorr(r+1:end,r+1:end); % only the passives

% correlation between replays
OL_OL = OL_corrs(triu(true(size(OL_corrs)),1));
PR_PR = PR_corrs(triu(true(size(PR_corrs)),1));

Cmat(pix,m+1) = mean(Cmat(pix,1:r)); % mean CL to OL
Cmat(pix,m+2) = mean(Cmat(pix,(r+1):m)); % mean CL to PR
Cmat(pix,m+3) = mean(OL_OL); % mean OL to OL
Cmat(pix,m+4) = mean(PR_PR); % mean PR to PR

end

%%
figure;
for i = 1:(r+p+4)
    subplot(4,5,i);
    imagesc(reshape(Cmat(:,i),stackdims(1),stackdims(2)));
    colormap(brewermap([50],'*RdBu'));
    set(gca,'XTick',[],'YTick',[]);
%     hold on
%     for x = 1:length(WhichROIs)
%         plot(WhichROIs(x,2),WhichROIs(x,1),'.k');
%     end
end

%% plot traces for selected ROIs across the 3 conditions
selROIs = [70 28; 70 38; 57 26; 57 45; 51 17; 35 25; 25 11; 31 11]; 
selROIs = vertcat(selROIs, fliplr(WhichROIs));
selROIs = sortrows(selROIs,[1 2]);
ts = frame_TS;
ts = ts(valid_frames)-ts(1);
%figure;
for roi = 1:length(selROIs)
    pix = stackdims(1)*(selROIs(roi,1)-1) + selROIs(roi,2);
    
    % get traces
    clear MyTraces nSamps
    CL_trace = WF_stack(pix,frame_idx(valid_frames));
    nSamps(1) = length(CL_trace);
    MyTraces(:,1) = CL_trace;
    for i = 1:r
        indices = Frames_idx(i,:);
        Replay_trace = WF_stack(pix,indices(1):indices(2));
        nSamps(i+1) = length(Replay_trace);
        MyTraces(1:length(Replay_trace),i+1) = Replay_trace;
    end
    for i = 1:p
        indices = WF_timestamps.Passive(Frames_idx(i+r,:),2);
        Replay_trace = WF_stack(pix,indices(1):indices(2));
        nSamps(i+r+1) = length(Replay_trace);
        MyTraces(1:length(Replay_trace),i+r+1) = Replay_trace;
    end
    
    MyTraces(1+min(nSamps):end,:) = [];

    meanOL  = mean(MyTraces(:,2:(r+1)),2);
    stdOL   = std(MyTraces(:,2:(r+1))');
    meanPR  = mean(MyTraces(:,(r+2):end),2);
    stdPR   = std(MyTraces(:,(r+2):end)');
    
    Lims = [10*floor(min(meanOL)/10) 10*ceil(max(meanOL)/10)];
    
    figure; 
    
    subplot(1,7,[2:7]);
    ProcessOpenLoopBasic(OpenLoop, SampleRate, TargetZones, 'plottrials', 1, 'TrialHeight',Lims);
    hold on
    plot(ts(1:min(nSamps)),MyTraces(:,1),'k');
    plot(ts(1:min(nSamps)),meanOL,'Color',Plot_Colors('t'));
    plot(ts(1:min(nSamps)),meanOL+stdOL',':','Color',Plot_Colors('t'));
    plot(ts(1:min(nSamps)),meanOL-stdOL',':','Color',Plot_Colors('t'));
    plot(ts(1:min(nSamps)),meanPR,'Color',Plot_Colors('r'));
    plot(ts(1:min(nSamps)),meanPR+stdPR',':','Color',Plot_Colors('r'));
    plot(ts(1:min(nSamps)),meanPR-stdPR',':','Color',Plot_Colors('r'));
    
    % lets also add here the single pixel correlation 
    
    refPix = WF_stack(pix,:); % entire timeseries
    clear C
    for j = 1:nPixels
        C(j,1) = corr(WF_stack(j,:)',refPix');
    end
    subplot(1,7,1);
    imagesc(reshape(C,stackdims(1),stackdims(2)));
    colormap(brewermap([50],'*RdBu'));
    hold on
    plot(selROIs(roi,1),selROIs(roi,2),'.k');
    title([num2str(selROIs(roi,1)),', ',num2str(selROIs(roi,2))]);
    set(gca,'XTick',[],'YTick',[]);
    
    set(gcf,'Position',[2055 331 1447  155]);
end

%% selected ROIs for the RO1
%% plot traces for selected ROIs across the 3 conditions
selROIs = [68 28; 51 17; 47 26; 37 26; 25 11]; 
ts = frame_TS;
ts = ts(valid_frames)-ts(1);
figure;
%clear C;
for roi = 1:length(selROIs)
    pix = stackdims(1)*(selROIs(roi,1)-1) + selROIs(roi,2);
    
    % get traces
    clear MyTraces nSamps
    CL_trace = WF_stack(pix,frame_idx(valid_frames));
    nSamps(1) = length(CL_trace);
    MyTraces(:,1) = CL_trace;
    for i = 1:r
        indices = Frames_idx(i,:);
        Replay_trace = WF_stack(pix,indices(1):indices(2));
        nSamps(i+1) = length(Replay_trace);
        MyTraces(1:length(Replay_trace),i+1) = Replay_trace;
    end
    for i = 1:p
        indices = WF_timestamps.Passive(Frames_idx(i+r,:),2);
        Replay_trace = WF_stack(pix,indices(1):indices(2));
        nSamps(i+r+1) = length(Replay_trace);
        MyTraces(1:length(Replay_trace),i+r+1) = Replay_trace;
    end
    
    MyTraces(1+min(nSamps):end,:) = [];

    meanOL  = mean(MyTraces(:,2:(r+1)),2);
    stdOL   = std(MyTraces(:,2:(r+1))');
    meanPR  = mean(MyTraces(:,(r+2):end),2);
    stdPR   = std(MyTraces(:,(r+2):end)');
    
    Lims = [10*floor(min(meanOL)/10) 10*ceil(max(meanOL)/10)];
        
    subplot(5,7,((roi-1)*7) + [2:7]);
    ProcessOpenLoopBasic(OpenLoop, SampleRate, TargetZones, 'plottrials', 1, 'TrialHeight',Lims);
    hold on
    plot(ts(1:min(nSamps)),MyTraces(:,1),'k');
    plot(ts(1:min(nSamps)),meanOL,'Color',Plot_Colors('t'));
%     plot(ts(1:min(nSamps)),meanOL+stdOL',':','Color',Plot_Colors('t'));
%     plot(ts(1:min(nSamps)),meanOL-stdOL',':','Color',Plot_Colors('t'));
    plot(ts(1:min(nSamps)),meanPR,'Color',Plot_Colors('r'));
%     plot(ts(1:min(nSamps)),meanPR+stdPR',':','Color',Plot_Colors('r'));
%     plot(ts(1:min(nSamps)),meanPR-stdPR',':','Color',Plot_Colors('r'));
    
    % lets also add here the single pixel correlation 
    
%     refPix = WF_stack(pix,:); % entire timeseries
%     for j = 1:nPixels
%         C(j,roi) = corr(WF_stack(j,:)',refPix');
%     end
    
    subplot(5,7,((roi-1)*7) + 1);
    imagesc(reshape(C(:,roi),stackdims(1),stackdims(2)),[-0.5 1]);
    colormap(brewermap([],'Greys'));
    hold on
    plot(selROIs(roi,1),selROIs(roi,2),'.k');
    title([num2str(selROIs(roi,1)),', ',num2str(selROIs(roi,2))]);
    set(gca,'XTick',[],'YTick',[]);
    
end
set(gcf,'Position',[2055 101 1447  5*155]);

%% just the mean image
WF_stack = reshape(WF_stack,stackdims(1),stackdims(2),nFrames);

