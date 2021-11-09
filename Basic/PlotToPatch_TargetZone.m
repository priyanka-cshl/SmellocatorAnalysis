function [handle_in] = PlotToPatch_TargetZone(handle_in, data_in, timestamp_in)
% data_in must be a column vector
a = data_in(:,1)'; b = data_in(:,2)';
% if any(a==500) || any(a==600)|| any(a==700)
%     b = ceil(b/100);
% end
if any(a>5)
    a(a>5) = 6;
    b(b~=0) = 5;
end
f = find(diff(a)~=0);
f = sort([1 length(a) f (f+1)]);
f = f';
X = [timestamp_in(f)' fliplr(timestamp_in(f)')];
Y = [a(f) fliplr(b(f))];
handle_in.Faces = 1:length(X);
handle_in.Vertices = [X' Y'];
end