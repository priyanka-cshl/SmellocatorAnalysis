
%FolderRoot = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/OB imaging/04-Aug-2022';
%FolderRoot = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/OB imaging/26-Jul-2022_4';
FolderRoot = '/mnt/data/OBImaging/testmouse3_brown/02-Sep-2022_2';
FolderRoot = '/mnt/data/OBImaging/testmouse2_black/07-Sep-2022';
[Ratio, TrialList] = GetRatioImages(FolderRoot);

%% load Trial List from the behavior side
if ~exist('TrialFile','var')
    TrialSequencePath = fullfile(FolderRoot,'OBimaging_20220907_o1.mat');
    %TrialSequencePath = fullfile(FolderRoot,'OBimaging_20220726_o0.mat');
else
    TrialSequencePath = fullfile(FolderRoot,TrialFile);
end

%%
load(TrialSequencePath);
TrialList_Behavior = session_data.TrialSequence;

% Assume first two trials in the Imaging side are weird
% check by fitting a line to motor locations from behavior side and imaging
% side (analog data)
ExtraTrials = size(TrialList,1) - size(TrialList_Behavior,1);

if ExtraTrials
    [~,gof] = fit(TrialList((ExtraTrials+1):end,2),TrialList_Behavior(:,1),'poly1');
    if gof.rsquare<0.9
        error('Unable to match Imaging and Behavior Trial Sequence');
    else
        Ratio(:,:,1:ExtraTrials) = [];
        TrialList(1:ExtraTrials,:) = [];
    end
end
Locations = unique(TrialList_Behavior(:,1));
Odors     = unique(TrialList_Behavior(:,2));
Reps      = size(TrialList_Behavior,1)/numel(Odors)/numel(Locations);

%%
for i = 1:numel(Odors)
    figure;
    colormap gray
    for j = 1:numel(Locations)
        % find all trials of a given odor and given location
        MyTrials = find(ismember(TrialList_Behavior,[Locations(j) Odors(i)],'rows'));
        
        for k = 1:numel(MyTrials)
            plotID = (k-1)*numel(Locations) + j;
            subplot(Reps,numel(Locations),plotID);
            imagesc(Ratio(:,:,MyTrials(k)));%,[0 0.1]);
            set(gca,'XTick',[],'YTick',[]);
        end
        
    end

end



    

