function [SniffsOut] = ParseSniffsByOccurence(AllSniffs, SortBy)

if nargin<2
    SortBy = 1;
end

%% split by stimulus conditions 
for snifftype = 1:3 % 3 odors

    % find the last sniff in the trial (odor was On and likely at center
    whichsniffs = intersect(find(AllSniffs(:,5)==snifftype), find(diff(AllSniffs(:,5))<0));

    switch snifftype
        case 1 % ITI
            whichsniffs = find(AllSniffs(:,4)==0);
        case 2
            whichsniffs = intersect(find(AllSniffs(:,4)==1), find(AllSniffs(:,5)==0));
        case 3
            whichsniffs = intersect(find(AllSniffs(:,4)==1), find(AllSniffs(:,5)==1));
        case 4
            whichsniffs = intersect(find(AllSniffs(:,4)==1), find(AllSniffs(:,5)==2));
        case 5
            whichsniffs = intersect(find(AllSniffs(:,4)==1), find(AllSniffs(:,5)==3));
    end

    SniffTS = AllSniffs(whichsniffs,:);   
    
    switch SortBy
        case 1 % by sessionphase and then duration
            SniffTS = sortrows(SniffTS,[8 3],'ascend');
        case 2
            % by sessionphase and then occurence
            SniffTS = sortrows(SniffTS,[8 1],'ascend');
    end
    
    SniffsOut{snifftype} = SniffTS;
end