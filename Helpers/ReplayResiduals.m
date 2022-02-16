function [ResidualsFull,Tags] = ReplayResiduals(PSTH,Trialcounts)

for MyUnit = 1:size(PSTH,3) % every unit
    MyFRs = squeeze(PSTH(:,:,MyUnit))'; % columns are different replays, rows are timepoints
    
    Residuals = []; Tags = [];
    % CL-OL
    whichCols = Trialcounts(1) + (1:Trialcounts(2)); % OL reps
    Residuals = vertcat(Residuals,mean((MyFRs(:,whichCols) - MyFRs(Trialcounts(1))).^2)');
    Tags      = vertcat(Tags,repmat(12,numel(whichCols),1));
    
    % CL-PR
    whichColsPassive = Trialcounts(1) + Trialcounts(2) + (1:Trialcounts(3));
    Residuals = vertcat(Residuals,mean((MyFRs(:,whichColsPassive) - MyFRs(Trialcounts(1))).^2)');
    Tags      = vertcat(Tags,repmat(13,numel(whichColsPassive),1));
    
    % OL-OL
    pairs = nchoosek(whichCols,2);
    Residuals = vertcat(Residuals, mean((MyFRs(:,pairs(:,1)) - MyFRs(:,pairs(:,2))).^2)');
    Tags      = vertcat(Tags,repmat(22,size(pairs,1),1));
    meanOL = mean(mean((MyFRs(:,pairs(:,1)) - MyFRs(:,pairs(:,2))).^2)');
    
    % PR-PR
    pairs = nchoosek(whichColsPassive,2);
    Residuals = vertcat(Residuals, mean((MyFRs(:,pairs(:,1)) - MyFRs(:,pairs(:,2))).^2)');
    Tags      = vertcat(Tags,repmat(33,size(pairs,1),1));
    
    % OL-PR
    pairs = [repmat(whichCols',numel(whichColsPassive),1) reshape(repmat(whichColsPassive,numel(whichCols),1),numel(whichCols)*numel(whichColsPassive),1)];
    Residuals = vertcat(Residuals, mean((MyFRs(:,pairs(:,1)) - MyFRs(:,pairs(:,2))).^2)');
    Tags      = vertcat(Tags,repmat(23,size(pairs,1),1));
    
    Residuals = Residuals/meanOL;
    
    ResidualsFull(:,MyUnit) = Residuals;
end