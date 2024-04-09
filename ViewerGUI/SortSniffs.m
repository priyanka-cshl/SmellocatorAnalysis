function [AllSniffs] = SortSniffs(AllSniffs,sortby)

switch sortby
    case 1
        % sort sniff List by Sniff Type, then sniff duration, then inh duration, then trial ID
        %AllSniffs = sortrows(AllSniffs,[4 15 16 1]);
        AllSniffs = sortrows(AllSniffs,[2 15 16 1]);
    case 2
        % sort sniff List by Sniff Type, then inh duration, then sniff duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[4 16 15 1]);
    case 3
        % sort by odor state, then odor location, then sniff duration, then inh duration, then trial ID
        AllSniffs = sortrows(AllSniffs,[2 13 15 16 1]);
    case 4
        % sort sniff List by Sniff Type, then prev-sniff duration, then current sniff duration, then trial ID
        % compute previous sniff duration
        % AllSniffs(:,17) = AllSniffs(:,7) - AllSniffs(:,5); this is already done in SelectSniffs.m
        AllSniffs = sortrows(AllSniffs,[4 17 15 1]);
end

end