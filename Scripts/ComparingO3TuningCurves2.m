%load '/mnt/data/DataMatrices/TuningCurves_O3_20211005.mat';


N = size(Curve_CL{1},4); % #units
% {} - different odors
% for each {}
% dimension 1 - bins
% dimension 2 - mean, median, sem, counts
% dimension 3 - actual, shuffled
% dimension 4 - units
Z = [];
for unit = 1:N
    for odor = 1:3
        X = [];
        a = max(Curve_CL{odor}(:,1,1,unit));
        X(:,1) = Curve_CL{odor}(:,1,1,unit)/a;
        X(:,2) = Curve_OL{odor}(:,1,1,unit)/a;
        X(:,3) = Curve_PR{odor}(:,1,1,unit)/a;
        X(:,4) = Curve_CL1{odor}(:,1,1,unit)/a;
        X(:,5) = Curve_CL2{odor}(:,1,1,unit)/a;
        
        %X([1:4 21],:) = []; % delete unsmapled locations
        X(1,:) = []; % delete unsmapled locations
        
        if any(isnan(X))
            keyboard;
        end
        
        % calculate residuals
        Y(1) = mean((X(:,1) - X(:,2)).^2); % CL-OL
        Y(2) = mean((X(:,1) - X(:,3)).^2); % CL-PR
        Y(3) = mean((X(:,2) - X(:,3)).^2); % OL-PR
        Y(4) = mean((X(:,4) - X(:,5)).^2); % CL-CL
                
        Z = vertcat(Z, [Y' (1:4)' unit*ones(4,1) odor*ones(4,1)]);
    end
end
