function [handle_in] = PlotToPatch_TrialTF(handle_in, trial_in, timestamp_in, offsets, TZ)

f = find(trial_in(:,1)>0);
data_in = zeros(length(trial_in),1);
data_in(f,1) = 1;
on_indices = timestamp_in( find(diff(data_in)==1) +1 );
off_indices = timestamp_in( find(diff(data_in)==-1) +1 );
while numel(off_indices)>numel(on_indices)
    off_indices(1,:) = [];
end
ValveTS = [on_indices off_indices]';
if ~isempty(ValveTS)
    handle_in.trial_tf.Vertices = [ ...
        reshape([ValveTS(:) ValveTS(:)]', 2*numel(ValveTS), []) , ...
        repmat([offsets(1) offsets(2) offsets(2) offsets(1)]',size(ValveTS,2),1)];
    handle_in.trial_tf.Faces = reshape(1:2*numel(ValveTS),4,size(ValveTS,2))';
    
    myTZ = [TZ(find(diff(data_in)==1) +2) TZ(find(diff(data_in)==1) +2)];
    VertexColor = [];
    for x = 1:size(myTZ,1)
        VertexColor(x,1) = -numel(abs(myTZ(x,1)):-0.0458:0);
        VertexColor(x,2) = numel(abs(myTZ(x,1)):0.0458:5);
        if myTZ(x,1)<0
            VertexColor(x,:) = -VertexColor(x,:);
        end
    end
    VertexColor(:,3:4) = VertexColor(:,[2 1]);
    VertexColor = VertexColor';
    handle_in.trial_tf.FaceVertexCData = VertexColor(:);
    handle_in.trial_tf.FaceColor = 'interp';
    %handle_in.trial_tf.FaceAlpha = 0.7;
end

end