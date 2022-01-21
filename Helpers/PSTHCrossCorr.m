function [C] = PSTHCrossCorr(PSTH,UnitIds)

catIdx = [];
for MyUnit = 1:size(PSTH,3) % every unit
    MyFRs = squeeze(PSTH(:,:,MyUnit))'; % columns are trials
    MyCorr = corrcoef(MyFRs);
    
    m = size(MyCorr,2) - 1;
    
    % correlation of replay to closed loop
    C(1:m,MyUnit) = MyCorr(1,2:end);
    
    MyCorr(1,:) = [];
    MyCorr(:,1) = [];
    
    % correlation between replays
    foo = MyCorr(triu(true(size(MyCorr)),1));
    n = length(foo) + m;
    
    C(m+1:n,MyUnit) = foo;
    
    catIdx = [catIdx; ones(m,1); zeros(n-m,1)];
end

UnitNames = cellfun(@num2str,num2cell(UnitIds(:)),'uniformoutput',false);

figure;
plotSpread(C,'categoryIdx',catIdx','categoryMarkers',{'*','o'},'categoryColors',{'r','k'},...
    'xNames', UnitNames);
legend('within replay','replay vs close loop')
