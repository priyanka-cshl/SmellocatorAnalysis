% function [] = GetSniffTuningSummary()
% end

MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q9/Q9_20221116_r0_processed.mat';
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/Q8/Q8_20221204_r0_processed.mat';

% load data
[handles.TrialAligned, handles.TrialInfo, ...
    handles.ReplayAligned, handles.ReplayInfo, ...
    handles.TuningAligned, handles.TuningInfo, ...
    handles.AllUnits] = ...
    PreprocessSpikesAndSniffs(MySession);

% get sniff time stamps and info for the sniffs we want to plot
[handles.SelectedSniffs] = SelectSniffs(handles.TrialAligned, handles.TrialInfo, [1 2 3], ...
                                'includeITI', 1);
                            
% sort the sniffs
for whichodor = 1:3
    handles.SelectedSniffs{whichodor} = ...
        SortSniffs(handles.SelectedSniffs{whichodor}, 1);
end

PSTHOffset = -1000;
for n = 1:size(handles.AllUnits.Spikes,2) % every unit
    for whichodor = 1:3
        [~,AllFR{n,whichodor}] = SniffAlignedPlot(handles.SelectedSniffs{whichodor}, handles.AllUnits.Spikes{n}, ...
            'plotevents', 0, ...
            'plotspikes', 0, ...
            'psth', 1, ...
            'psthoffset', PSTHOffset, ...
            'alignto', 1, ...
            'warptype', 0);
    end
end

%%
for n = 1:size(handles.AllUnits.Spikes,2) % every unit
    
    for whichodor = 1:3
        FR = AllFR{n,whichodor};
        t = 4;
        myFR = (FR{t} - mean(FR{t}(451:500,1)))/mean(FR{t}(451:500,1));
        myFR(1:500,:) = [];
        l = numel(myFR) - 50;
        Tuning{whichodor}(n,1:numel(myFR)) = myFR;
        [p, idx] = max((myFR(1:l)));
        Peak(n,1:2,whichodor) = [myFR(idx) idx];
    end
    
end

%%
figure
for whichodor = 1:3
    Peak2 = Peak(:,:,whichodor);
    Peak2(:,3) = Peak2(:,1)>=0;
    Peak2(:,4) = 1:size(Peak2,1);
    Peak2 = sortrows(Peak2,[3 2]);
    subplot(1,3,whichodor)
    imagesc(Tuning{whichodor}(Peak2(:,4),:),[-1 5]);
    colormap(brewermap([100],'*RdBu'));
    set(gca,'TickDir','out');
end
