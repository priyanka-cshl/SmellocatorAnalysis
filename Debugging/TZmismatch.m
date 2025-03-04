% debugging the mismatch between suggested and actual Targtzone limits

FileName = '/mnt/grid-hs/pgupta/Behavior/Q5/Q5_20221209_r0.mat';

% load the data
Temp = load(FileName,'session_data');
MyTraces = Temp.session_data.trace; % data

trialVec = MyTraces(:,7);
trialVec(trialVec>0) = 1;
TrialOns = find(diff(trialVec)==1);

buggytrial = 105; % 104 from Ryan
% get data indices for this, prev and next trial
for t  = -1:1:1
    trialIdx = TrialOns((buggytrial+t))+1;
    trialIdx(2) = find(trialVec(trialIdx:end)==0,1,'first') + trialIdx - 2;
    idxs{t+2} = trialIdx(1):trialIdx(2);
end

figure;
for t = 1:3
    subplot(4,3,t); % lever DAC vs. raw
    hold on
    plot(MyTraces(idxs{t},1),MyTraces(idxs{t},2),'ok');
    set(gca,'XLim',[0 5], 'YLim',[0 3])

    subplot(4,3,t+3); % RE vs. inferred position
    hold on
    plot(MyTraces(idxs{t},3),MyTraces(idxs{t},4),'or');
    set(gca,'XLim',[0 4], 'YLim',[-100 100])

    subplot(4,3,t+6); % lever vs. RE
    hold on
    plot(MyTraces(idxs{t},1),MyTraces(idxs{t},3),'ob');
    set(gca,'XLim',[0 5], 'YLim',[0 4])

    subplot(4,3,t+9); % lever vs. RE
    hold on
    plot(MyTraces(idxs{t},1),MyTraces(idxs{t},4),'og');
    set(gca,'XLim',[0 5], 'YLim',[-100 100])

end


