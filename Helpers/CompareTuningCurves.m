function [AllCorrs,AllResiduals,AllControlCorrs,AllControlResiduals] = CompareTuningCurves(TuningCurves)
% dims: bins, mean sem, unit, odor

AllCorrs = []; AllResiduals = []; AllControlCorrs = []; AllControlResiduals = [];
for i = 1:size(TuningCurves.ClosedLoopFull,3) % every unit
    % stack CL, OL, PR together
    thisUnitCurves = squeeze([TuningCurves.ClosedLoopFull(:,1,i,:) ...
                              TuningCurves.OpenLoop(:,1,i,:) ...
                              TuningCurves.Passive(:,1,i,:)]);
    
    thisUnitControlCurves = squeeze([TuningCurves.ClosedLoopHalf1(:,1,i,:) ...
                                     TuningCurves.ClosedLoopHalf2(:,1,i,:)]);                      
    
    for odor = 1:size(thisUnitCurves,3) % every odor
        thisOdorCurves = thisUnitCurves(:,:,odor);
        % delete NaNs
        thisOdorCurves(find(isnan(sum(thisOdorCurves,2))),:) = [];
        % normalize to the max of this unit's tuning curves
        thisOdorCurves = thisOdorCurves/max(thisUnitCurves(:));
        
        % calculate correlations
        foo = corrcoef(thisOdorCurves);
        thisOdorCorrs = [foo(2,1) foo(3,1) foo(3,2)]; % CL-OL CL-PR OL-PR
        thisOdorResiduals = [mean(diff(thisOdorCurves(:,[1 2])')'.^2) ...
                             mean(diff(thisOdorCurves(:,[1 3])')'.^2) ...
                             mean(diff(thisOdorCurves(:,[2 3])')'.^2)];
        
                         
        AllCorrs = vertcat(AllCorrs, [thisOdorCorrs i odor]);
        AllResiduals = vertcat(AllResiduals, [thisOdorResiduals i odor]);
        
        thisOdorControlCurves = thisUnitControlCurves(:,:,odor);
        % delete NaNs
        thisOdorControlCurves(find(isnan(sum(thisOdorControlCurves,2))),:) = [];
        % normalize to the max of this unit's tuning curves
        thisOdorControlCurves = thisOdorControlCurves/max(thisUnitCurves(:));
        
        % calculate correlations
        foo = corrcoef(thisOdorControlCurves);
        thisOdorControlCorrs = foo(2,1); % CL-OL CL-PR OL-PR
        thisOdorControlResiduals = mean(diff(thisOdorControlCurves(:,[1 2])')'.^2);
                         
        AllControlCorrs = vertcat(AllControlCorrs, [thisOdorControlCorrs i odor]);
        AllControlResiduals = vertcat(AllControlResiduals, [thisOdorControlResiduals i odor]);
        
    end
    
end
end