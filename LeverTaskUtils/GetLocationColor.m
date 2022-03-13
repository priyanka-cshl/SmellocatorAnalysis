function [whichcolor] = GetLocationColor(x)
CMAP = colormap(brewermap([100],'rdbu'));
%set(gca,'CLim',120*[-1 1]);
foo = linspace(-120, 120, 100);
[y,z] = min(abs(foo-x));
whichcolor = CMAP(z(1),:);
end