function [handle_in] = PlotToPatch_Trial(handle_in, trial_in, timestamp_in, offsets, whichcase)
if nargin<5
    whichcase = 0;
end
% data_in must be a column vector
for i = 1:4 % three odors
    f = find(trial_in(:,1)==i);
    data_in = zeros(length(trial_in),1);
    data_in(f,1) = 1;
    on_indices = timestamp_in( find(diff(data_in)==1) +1 );
    off_indices = timestamp_in( find(diff(data_in)==-1) +1 );
    while numel(off_indices)>numel(on_indices)
        off_indices(1,:) = [];
    end
    ValveTS = [on_indices off_indices]';
    if ~isempty(ValveTS)
        handle_in.(['trial_on_',num2str(i)]).Vertices = [ ...
            reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
            repmat([offsets(1) offsets(2) offsets(2) offsets(1)]',size(ValveTS,2),1)];
        handle_in.(['trial_on_',num2str(i)]).Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
    end
end
end