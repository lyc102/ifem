function [elem,idx] = fixorientation(node,elem)
%% FIXORIENTATION set all triangles counter-clockwise
% 
%   elem = fixorientation(node,elem) computes signed area of all triangles
%   in the triangulation and switch the vertices such that all signed area
%   is positive, i.e. the orientation of all triangles are
%   counter-clockwise.
%
%
% See also fixorder, fixorder3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = -ve3(:,1).*ve2(:,2)+ve3(:,2).*ve2(:,1);
idx = find(area<0); 
elem(idx,[2 3]) = elem(idx,[3 2]);