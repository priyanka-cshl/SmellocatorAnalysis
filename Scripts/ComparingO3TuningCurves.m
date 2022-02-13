%load '/mnt/data/DataMatrices/TuningCurves_O3_20211005.mat';


N = size(Curve_CL{1},4); % #units
% {} - different odors
% for each {}
% dimension 1 - bins
% dimension 2 - mean, median, sem, counts
% dimension 3 - actual, shuffled
% dimension 4 - units

counts = 0;
for unit = 1:N
    if mod(unit,24) == 1
        figure;
    end
    XX = [];
    for odor = 1:3
        X = [];
        a = max(Curve_CL{odor}(:,1,1,unit));
        X(:,1) = Curve_CL{odor}(:,1,1,unit)/a;
        X(:,2) = Curve_OL{odor}(:,1,1,unit)/a;
        X(:,3) = Curve_PR{odor}(:,1,1,unit)/a;
        
        %X([1:4 21],:) = []; % delete unsmapled locations
        X(1,:) = []; % delete unsmapled locations
        
        if any(isnan(X))
            keyboard;
        end
        
        % calculate correlations
        foo = corrcoef(X);
        Y(:,1,unit,odor) = [foo(2,1) foo(3,1) foo(3,2)]; % CL-OL, CL-PR, OL-PR
        
        % calculate residuals
        Y(1,2,unit,odor) = mean((X(:,1) - X(:,2)).^2); % CL-OL
        Y(2,2,unit,odor) = mean((X(:,1) - X(:,3)).^2); % CL-PR
        Y(3,2,unit,odor) = mean((X(:,2) - X(:,3)).^2); % OL-PR
        % X(:,3,unit,odor)
        XX = vertcat(XX,[NaN NaN NaN],...
            [Curve_CL{odor}(:,1,1,unit) Curve_OL{odor}(:,1,1,unit) Curve_PR{odor}(:,1,1,unit)]);
        
        
        Z(counts+1:counts+3,:) = [Y(:,1,unit,odor) Y(:,2,unit,odor) (1:3)' unit*ones(3,1) odor*ones(3,1)];
        counts = counts + 3;
    end
    % lets plot some stuff
    whichplot = mod(unit,24);
    if ~whichplot
        whichplot = 24;
    end
    subplot(4,6,whichplot); hold on
    plot(XX);
    title([num2str(foo(2,1)),',',num2str(foo(3,1)),',',num2str(foo(3,2))])
end


figure;
for odor = 1:3
    for conditions = 1:3
        subplot(3,3,odor+3*(conditions-1));
        xvar = Y(conditions,1,:,odor); % correlation
        yvar = Y(conditions,2,:,odor); % residuals
        scatter(xvar,yvar);
    end
end