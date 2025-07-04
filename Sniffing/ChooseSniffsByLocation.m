function [SniffsOut] = ChooseSniffsByLocation(AllSniffs, SortBy, timechunks)

if nargin<2
    SortBy = 1;
    timechunks = [];
elseif nargin<3
    timechunks = [];
end

targetLocs = 8;

%% split by stimulus conditions 
for snifftype = 1:5

    switch snifftype
        case 1 % ITI
            %whichsniffs = find(AllSniffs(:,4)==0);
            whichsniffs = intersect(find(AllSniffs(:,4)==0), find(AllSniffs(:,5)>=0));
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

    if snifftype<3 % not odors
        % choose the non-target sniffs
        whichsniffs = find(SniffTS(:,6)<=-2*targetLocs | SniffTS(:,6)>=2*targetLocs);
    else
        % choose the target sniffs
        whichsniffs = find(SniffTS(:,6)>=-targetLocs & SniffTS(:,6)<=targetLocs);
    end
    SniffTS = SniffTS(whichsniffs,:);
    
    switch SortBy
        case 1 % by sessionphase and then duration
            SniffTS = sortrows(SniffTS,[8 3],'ascend');
        case 2
            % by sessionphase and then occurence
            SniffTS = sortrows(SniffTS,[8 1],'ascend');
        case 3
            % by occurence
            SniffTS = sortrows(SniffTS,[1],'ascend');

            % and label by recording chunks
            for chunks = 1:size(timechunks,1)
                SniffTS(find(SniffTS(:,1)>timechunks(chunks,1)),8) = chunks;
            end
    end
    
    SniffsOut{snifftype} = SniffTS;
end