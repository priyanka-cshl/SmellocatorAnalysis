function [handle_in] = PlotToPatch_Perturbation(handle_in, trial_in, timestamp_in, offsets)

f = find(trial_in(:,1)~=0);
data_in = zeros(length(trial_in),1);
data_in(f,1) = 1;
on_indices = timestamp_in( find(diff(data_in)==1) +1 );
off_indices = timestamp_in( find(diff(data_in)==-1) +1 );
perturbValue = trial_in( find(diff(data_in)==1) +1 );
ValveTS = [on_indices off_indices]';
if ~isempty(ValveTS)
    handle_in.perturb.Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat([offsets(1) offsets(2) offsets(2) offsets(1)]',size(ValveTS,2),1)];
    handle_in.perturb.Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
    
    VertexColor = [];
    for x = 1:numel(perturbValue)
        VertexColor(x,1) = perturbValue(x);
        VertexColor(x,2) = perturbValue(x);
    end
    VertexColor(:,3:4) = VertexColor(:,[2 1]);
    VertexColor = VertexColor';
    handle_in.perturb.FaceVertexCData = VertexColor(:);
    handle_in.perturb.FaceColor = 'interp';
    %handle_in.perturb.FaceAlpha = 0.7;
end

end

% test 