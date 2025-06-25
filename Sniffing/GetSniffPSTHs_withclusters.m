function [] = GetSniffPSTHs_withClusters(procFile)
%% get all PSTHs for all units
%  get correlations within and across units
%  use the .wdw processed files

baseFolder = '/home/priyanka/Desktop/forWDW';
%procFile = 'Q9_20221119_r0_processed.mat';
% procFile = 'Q4_20221112_r0_processed.mat';
% procFile = 'Q8_20221209_r0_processed.mat';
% procFile = 'Q5_20221122_r0_processed.mat';
% procFile = 'NS12_20230727_r0_processed.mat';
procFile = 'nO3_20211005_r0_processed.mat';
%procFile = 'Q5_20221122_r0_processed.mat';
%load(fullfile(baseFolder,procFile));

% get a list of sniffs, also get units
[AllSniffs, ~, SingleUnits] = GetAllSniffs(fullfile(baseFolder,procFile));

% get sniff time stamps and info for the sniffs we want to plot
[ParsedSniffs] = ChooseSniffsByLocation(AllSniffs,2);

%% get PSTHs
nUnits = size(SingleUnits,2);
% psth settings
psthbins = 20;
sniffwindow = [-0.1 0.5];
PSTHbig = []; 
for i = 1:nUnits
    thisUnitspikes = SingleUnits(i).spikes;
    [~, SniffPSTHs, ~] = GetSniffLockedSpikesWithPSTH(ParsedSniffs, thisUnitspikes, 'PSTHBinsize', psthbins, 'window', sniffwindow, 'onlyCL', 1);
    PSTHbig(:,:,i) = (cell2mat(cellfun(@(x)(x.mean), SniffPSTHs, 'UniformOutput', false)'))';
end

% PSTHraw = reshape(PSTHbig,size(PSTHbig,1),size(PSTHbig,2)*size(PSTHbig,3));
% % subtract ITI 
% PSTHsub = reshape(PSTHbig-PSTHbig(:,1,:),size(PSTHbig,1),size(PSTHbig,2)*size(PSTHbig,3));
%PSTHsmooth = sgolayfilt(PSTHbig,1,3);   
end
