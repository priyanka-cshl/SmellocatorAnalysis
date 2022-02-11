function [Residuals] = ReplayResiduals(PSTH,Trialcounts)

for MyUnit = 1:size(PSTH,3) % every unit
    MyFRs = squeeze(PSTH(:,:,MyUnit))'; % columns are different replays, rows are timepoints

    % within open loop
    Residual = []; ResidualPassive = []; Residual_CL = []; Residual_PR = []; Residual_OL_PR = [];
    counts = 0;
    whichCols = Trialcounts(1) + (1:Trialcounts(2));
    whichColsPassive = Trialcounts(1) + Trialcounts(2) + (1:Trialcounts(3));
    for i = 1:numel(whichCols)
        x = whichCols(i);
        y = whichCols(whichCols~=x);
        % take the mean of all but this Col and calculate residual w.r.t. to it
        Residual(i) = mean((mean(MyFRs(:,y),2) - MyFRs(:,x)).^2); %#ok<AGROW>
        % also calculate residual w.r.t. to closed-loop
        Residual_CL(i) = mean((mean(MyFRs(:,y),2) - MyFRs(:,1)).^2); %#ok<AGROW>

        % also calculate residual w.r.t. Passive
        for j = 1:numel(whichColsPassive)
            if i == 1
                a = whichColsPassive(j);
                b = whichColsPassive(whichColsPassive~=a);
                % take the mean of all but this Col and calculate residual w.r.t. to it
                ResidualPassive(j) = mean((mean(MyFRs(:,b),2) - MyFRs(:,a)).^2); %#ok<AGROW>
                % also calculate residual w.r.t. to closed-loop
                Residual_PR(j) = mean((mean(MyFRs(:,b),2) - MyFRs(:,1)).^2); %#ok<AGROW>
            end
            % residual of a given passive trial to openloop
            counts = counts + 1;
            Residual_OL_PR(counts) = mean((mean(MyFRs(:,y),2) - MyFRs(:,a)).^2); %#ok<AGROW>
        end
    end
    Residuals(MyUnit).withinOL = [mean(Residual) std(Residual)]; %#ok<*AGROW>
    Residuals(MyUnit).CL2OL = [mean(Residual_CL) std(Residual_CL)];
    Residuals(MyUnit).withinPR = [mean(ResidualPassive) std(ResidualPassive)];
    Residuals(MyUnit).CL2PR = [mean(Residual_PR) std(Residual_PR)];
    Residuals(MyUnit).PR2OL = [mean(Residual_OL_PR) std(Residual_OL_PR)];
    
end