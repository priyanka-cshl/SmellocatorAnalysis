
%load(fullfile(myKsDir,'ClustersFull.mat')) %,"cluster",'-v7.3')

nUnits = size(cluster,2);
rounds = size(cluster(1).WF,2);
nPSTH = size(cluster(1).SniffPSTH,2);
WFLength = size(cluster(1).WF{1},3);
AllWFs = [];
AllPSTHs = [];

for whichunit = 1:nUnits
    tetrodeChannels = (floor(cluster(whichunit).tetrode)-1)*4 + (1:4);
    tetrodeChannels = 1:40;
    concatWF = []; 
    for i = 1:rounds
        concatWF(:,i) = reshape(...
            squeeze(cluster(whichunit).meanWF{i}(tetrodeChannels,:,1))' ,...
            WFLength*numel(tetrodeChannels), 1);
    end

    AllWFs = horzcat(AllWFs, concatWF);
    
    concatPSTH = [];
    for j = 1:nPSTH
        concatPSTH = horzcat(concatPSTH, ...
            cluster(whichunit).SniffPSTH{1,j}.mean(:,1:end-1) );
    end
    
    AllPSTHs = horzcat(AllPSTHs, concatPSTH');

end

WFCrossCorr = corr(AllWFs);
SniffCrossCorr = corr(AllPSTHs);

