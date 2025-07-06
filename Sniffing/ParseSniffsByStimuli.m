function [SniffsOut,odors] = ParseSniffsByStimuli(AllSniffs, varargin)

%% extract inputs
narginchk(1,inf);
params = inputParser;
params.CaseSensitive = false;
params.addParameter('SortBy', 1, @(x) isnumeric(x));
params.addParameter('SelectLocs', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
SortBy = params.Results.SortBy;
SelectLocs = params.Results.SelectLocs;

odors = unique(AllSniffs(:,5)); % different odors

%% split by stimulus conditions 
for snifftype = 1:(numel(odors)+1)

    switch snifftype
        case 1 % ITI
            whichsniffs = find(AllSniffs(:,4)==0);
        otherwise
            whichsniffs = intersect(find(AllSniffs(:,4)==1), find(AllSniffs(:,5)==odors(snifftype-1)));
    end

    SniffTS = AllSniffs(whichsniffs,:);   
    
    switch SortBy
        case 1 % by sessionphase and then duration
            SniffTS = sortrows(SniffTS,[8 3],'ascend');
        case 2
            % by sessionphase and then occurence
            SniffTS = sortrows(SniffTS,[8 1],'ascend');
        case 3
            % by sessionphase and then location, then duration
            SniffTS = sortrows(SniffTS,[8 6 3],'ascend');
    end
    
    SniffsOut{snifftype} = SniffTS;
end