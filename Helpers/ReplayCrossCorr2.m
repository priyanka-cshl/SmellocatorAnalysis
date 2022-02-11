function [C,g] = ReplayCrossCorr2(PSTH,Trialcounts)

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

    C(MyUnit).withinOL = [median(OL_OL) std(OL_OL)];
    C(MyUnit).CL2OL = [median(CL_OL) std(CL_OL)];
    C(MyUnit).withinPR = [median(PR_PR) std(PR_PR)];
    C(MyUnit).CL2PR = [median(CL_PR) std(CL_PR)];
    C(MyUnit).PR2OL = [median(OL_PR) std(OL_PR)];
end
