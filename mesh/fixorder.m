function [elem,idx,area] = fixorder(node,elem)
%% FIXORDER set all triangles counter-clockwise
% 
%   elem = FIXORDER(node,elem) computes signed area of all triangles in the
%   triangulation and switch the vertices such that all signed area is
%   positive, i.e. the orientation of all triangles are counter-clockwise.
%
%   [elem,idx,area] = FIXORDER(node,elem) also outputs the index set of
%   elements whose area is negative and the absolute value of area.
%
%   fixorder is recommend to use if the mesh is obtained by delaunay
%   function in matlab.
%
% See also fixorientation, fixorder3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

% compute signed area of each triangle
area = simplexvolume(node,elem);
% find triangles with negative area and switch the vertices
idx = find(area<0); 
elem(idx,[2 3]) = elem(idx,[3 2]);
area(idx) = -area(idx);