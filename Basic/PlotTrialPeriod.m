function [] = PlotTrialPeriod(TrialDeets)

%TrialDeets = n x 3; n = no. of trials, each row = on and off indices, TZ value
TS = TrialDeets(:,1:2)';
myTZ = TrialDeets(:,3);
offsets = [0 5];

myVertices  = [ reshape([TS(:) TS(:)]', 2*numel(TS), []) , ...
        repmat([offsets(1) offsets(2) offsets(2) offsets(1)]',size(TS,2),1)] ;

myFaces     = reshape(1:2*numel(TS),4,size(TS,2))';
myVertexColor = [];
for x = 1:size(TrialDeets,1)
    myVertexColor(x,[1 4]) = -numel(abs(myTZ(x,1)):-0.0458:0);
    myVertexColor(x,[2 3]) = numel(abs(myTZ(x,1)):0.0458:5);
end
myVertexColor = myVertexColor';

fill(myVertices(:,1),myVertices(:,2),[1 1 0],...
    'EdgeColor','none',...
    'FaceColor','interp','Faces',myFaces,'FaceVertexCdata',myVertexColor(:));

colormap(brewermap([100],'rdbu'));
set(gca,'CLim',120*[-1 1]);

end