
%% function to plot trials and segmented move analysis in one figure
BasePath = '/home/priyanka/Dropbox';
BasePath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator';

for m = 1; %1:5
    switch m
        case 1
            WhereSession = 'Smellocator_Behavior_2/Q3/Q3_20221012_r0_processed.mat';
        case 2
            WhereSession = 'Smellocator_Behavior_2/Q4/Q4_20221129_r0_processed.mat';
        case 3
            WhereSession = 'Smellocator_Behavior_2/Q5/Q5_20221028_r0_processed.mat';
        case 4
            WhereSession = 'Smellocator_Behavior_2/Q8/Q8_20221122_r0_processed.mat';
        case 5
            WhereSession = 'Smellocator_Behavior_2/Q9/Q9_20221108_r0_processed.mat';
    end
    WhereSession = fullfile(BasePath,WhereSession);

    clearvars -except m BasePath WhereSession fitstats xystats
    %%
    % load the data
    [TracesOut, whichTraces, SegmentedLever, thisSniffParams, TrialTimeStamps, TrialIndices, startspositive, ...
        TrialInfo, SampleRate, TargetZones] = ...
        LoadProcessedLeverSession(WhereSession);
    [~,MouseName] = fileparts(fileparts(WhereSession));

    %% per trial
    for whichtrial = 1:size(TrialIndices,1)
        % segments for this trial
        whichsniffs = find(thisSniffParams(:,11)==mean(TrialTimeStamps(whichtrial,:)));
        x_vals = thisSniffParams(whichsniffs,1);
        if ~startspositive
            x_vals = -x_vals;
        end
        y_vals = thisSniffParams(whichsniffs,4);
        
        if numel(x_vals)>=2
            % fit a line?
            [f,gof] = fit(x_vals,y_vals,'poly1');
            rmse = [x_vals*f.p1 + f.p2 - y_vals].^2;
            fitstats.(MouseName)(whichtrial,:) = [f.p1 f.p2 gof.rsquare gof.adjrsquare gof.rmse mean(rmse) median(rmse) (median(rmse) - mean(rmse))];
            xystats.(MouseName)(whichtrial,:) = [median(x_vals) std(x_vals) median(y_vals) std(y_vals)];
        end

    end

    %% plot trials sorted by some metric?
    % colors
    OdorColors(1,:) = [.8 .8 .8];
    OdorColors(2,:) = [0.8941    0.9412    0.9020];
    OdorColors(3,:) = [0.8706    0.9216    0.9804];
    OdorColors(4,:) = [0.93    0.84    0.84];
    TZColor(1,:)    = [1 0 0];
    TZColor(2,:)    = [1 1 0];
    LeverInh(1,:)   = [199 21 133]./256;
    handles.plotcolors.resp     = [52 101 164]./256;
    handles.plotcolors.licks    = [199 21 133]./256; %[239 41 41]./256;
    handles.plotcolors.red      = [239 41 41]./256;
    handles.plotcolors.rewards  = [0 139 139]./256;
    handles.plotcolors.LeverSeg = [173 127 168]./256; %[0 0 0];

    figure; 
    set(gcf, 'color', 'w');
    trial_gap = 0.5;
    time_offset = trial_gap;
    trialmax = 20;
    nplots = 11;
    p1 = 1;
    plotmax = 100;

    pooledstats = [(1:size(fitstats.(MouseName),1))'  fitstats.(MouseName) xystats.(MouseName)];
    % trial# slope p2 r2 adjr2 rmse mean-rmse med-rmse diff-rmse (1:9)
    % median-loc std-loc median-disp std-disp (10:13)
    %[~,Trialorder] = sortrows([fitstats.(MouseName) xystats.(MouseName)],[1 9]); % by slope
    %pooledstats = sortrows(pooledstats,[7 2]); % slope, diff rmse
    pooledstats = sortrows(pooledstats,[7]); % slope, diff rmse
    Trialorder = pooledstats(:,1);
    trial_durations = TrialTimeStamps(:,2)-TrialTimeStamps(:,1);
    trial_durations = trial_durations(Trialorder);

    max_duration = 0; 
    for i = 1:trialmax:numel(trial_durations)
        j = min(i+(trialmax-1),numel(trial_durations));
        max_duration = max(max_duration,sum(trial_durations(i:j)));
    end
    max_duration = max_duration + trialmax*trial_gap;

    % figure out trial chunks
    foo = trial_gap + cumsum(trial_durations) + cumsum(repmat(trial_gap,numel(trial_durations),1));
    foo(end+1) = foo(end) + trial_gap;
    nextlap = 130; lapsize = nextlap;
    for thistrial = 1:size(TrialIndices,1)
        %whichtrial = Trialorder(thistrial);
        whichtrial = pooledstats(thistrial,1);

        ts = TrialTimeStamps(whichtrial,:);
        ts = ts - ts(1) + time_offset;

%         if mod(thistrial,trialmax) == 1
%             xlimstart = ts(1);
%             p1 = p1+1;
%             subplot(nplots,1,p1); hold on
%         end
        if thistrial == 1
            xlimstart = ts(1);
            subplot(nplots,1,p1); hold on
        end

        if foo(thistrial)>nextlap  
            xlimstart = ts(1);
            p1 = p1+1;
            nextlap = nextlap + lapsize;
            subplot(nplots,1,p1); hold on
        end

        subplot(nplots,1,p1);
        % plot the odor box
        rectangle('Position',[ts(1) 0 diff(ts) 5],'EdgeColor','none','FaceColor',OdorColors(TrialInfo.Odor(whichtrial),:));
        % plot the target zone
        tzlims = TargetZones(TrialInfo.TargetZoneType(whichtrial),[3 1]);
        patch('Faces',1:4,'Vertices',[ [ts(1) ts(1) ts(2) ts(2)]' [tzlims(1) tzlims(2) tzlims(2) tzlims(1)]'], ...
            'EdgeColor','none','FaceColor',TZColor(TrialInfo.Success(whichtrial)+1,:),'FaceAlpha',0.2);
        % plot lever traces
        idx = TrialIndices(whichtrial,:); % trial start and end
        idx(1) = idx(1) - round(SampleRate*trial_gap/2);
        timevals = TracesOut(idx(1):idx(2),4) - TrialTimeStamps(whichtrial,1) + time_offset;
        plot(timevals, TracesOut(idx(1):idx(2),1), 'k', 'LineWidth', 2);
        plot(timevals, TracesOut(idx(1):idx(2),12), 'color', LeverInh, 'LineWidth', 2);
        % trial number
        text(mean(ts),4.8,num2str(whichtrial),'FontSize',6,'HorizontalAlignment','center');

        % keep track of offset
        time_offset = ts(2) + trial_gap;

%         if mod(thistrial,trialmax) == 0
%             subplot(nplots,1,p1);
%             set(gca,'XTick',[],'YTick',[], 'XLim', xlimstart + [0 max_duration]);
% 
%             if mod(p1,nplots) == 0
%                 set(gcf,'Position',[2000 50 2000 1000]);
%                 %figname = [SessionName,'_',num2str(whichtrial)];
%                 %saveas(gcf,fullfile(FigPath,[figname,'.png']));
%                 %close(gcf);
%                 p1 = 0;
%                 figure;
%             end
% 
%         end
        if thistrial == size(TrialIndices,1)
            subplot(nplots,1,p1);
            set(gca,'XTick',[],'YTick',[], 'XLim', xlimstart + [-trial_gap lapsize+(2*trial_gap)]);
            set(gcf,'Position',[2000 50 2000 1000]);

        elseif foo(thistrial+1)>nextlap
            subplot(nplots,1,p1);
            set(gca,'XTick',[],'YTick',[], 'XLim', xlimstart + [-trial_gap lapsize+(2*trial_gap)]);

            if mod(p1,nplots) == 0
                set(gcf,'Position',[2000 50 2000 1000]);
                %figname = [SessionName,'_',num2str(whichtrial)];
                %saveas(gcf,fullfile(FigPath,[figname,'.png']));
                %close(gcf);
                p1 = 1;
                figure;
                set(gcf, 'color', 'w');
            end

        end
    end
end

