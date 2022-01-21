function [C,g] = ReplayCrossCorr(PSTH,Trialcounts)

for MyUnit = 1:size(PSTH,3) % every unit
    MyFRs = squeeze(PSTH(:,:,MyUnit))'; % columns are trials
    FullCorr = corrcoef(MyFRs);
    
    % correlation of active replay to closed loop
    CL_OL = FullCorr(1,1+(Trialcounts(1):Trialcounts(2)));
    
    % correlation of passive replay to closed loop
    CL_PR = FullCorr(1,(1+Trialcounts(1)+Trialcounts(2)):end);
    
    % within OL correlation
    MyCorr = FullCorr(1+(Trialcounts(1):Trialcounts(2)),1+(Trialcounts(1):Trialcounts(2)));
    OL_OL = MyCorr(triu(true(size(MyCorr)),1));
    
    % within PR correlation
    MyCorr = FullCorr((1+Trialcounts(1)+Trialcounts(2)):end,(1+Trialcounts(1)+Trialcounts(2)):end);
    PR_PR = MyCorr(triu(true(size(MyCorr)),1));
    
    % OL to PR correlation
    MyCorr = FullCorr(1+(Trialcounts(1):Trialcounts(2)),(1+Trialcounts(1)+Trialcounts(2)):end);
    OL_PR = MyCorr(:);
    
    C(:,MyUnit) = vertcat(CL_OL',CL_PR',OL_OL,PR_PR,OL_PR);
end

g1 = repmat({'CL-OL'},numel(CL_OL),1);
g2 = repmat({'CL-PR'},numel(CL_PR),1);
g3 = repmat({'OL-OL'},numel(OL_OL),1);
g4 = repmat({'PR-PR'},numel(PR_PR),1);
g5 = repmat({'OL-PR'},numel(OL_PR),1);
g = [g1; g2; g3; g4; g5];
