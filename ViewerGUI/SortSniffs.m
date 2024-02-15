function [AllSniffs] = SortSniffs(AllSniffs,sortby)

switch sortby
    case 1
        % sort sniff List by Sniff Type, then sniff duration, then inh duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[4 15 16 1]);
    case 2
        % sort sniff List by Sniff Type, then inh duration, then sniff duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[4 16 15 1]);
    case 3
        % sort by odor state, then odor location, then sniff duration, then inh duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[2 13 15 16 1]);
end

end