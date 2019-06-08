function showmeshdensity(node,elem)
%% SHOWMESHDENSITY show the density of the mesh
%
% showmeshdensity(node,elem) display the size of each triangle using different
% color.
%
% Example
%   load lakemesh
%   showmeshdensity(node,elem)
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*abs(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
h = sqrt(area);
showsolution(node,elem,h)
colorbar
axis equal; axis tight; axis off;