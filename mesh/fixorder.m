function [elem,idx,area,bdFlag] = fixorder(node,elem,bdFlag)
%% FIXORDER set all triangles counter-clockwise
% 
%   elem = FIXORDER(node,elem) computes signed area (volume) of all
%   triangles (tetrahedron) in the triangulation and switch the vertices
%   such that all signed area (volume) are positive.
%
%   [elem,idx,area] = FIXORDER(node,elem) also outputs the index set of
%   elements whose area (volume) are negative and the absolute value of area.
%
%   [elem,idx,area,bdFlag] = FIXORDER(node,elem,bdFlag) changes the bdFlag
%   for boundary conditions. 
%
%   fixorder is recommend to use if the mesh is obtained by delaunay
%   function in matlab.
%
% See also fixorientation, fixorder3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin==2)
    bdFlag = [];
end
% compute signed area of each triangle
[area,elemSign] = simplexvolume(node,elem);
% find triangles with negative area and switch the vertices
idx = find(elemSign==-1);
elem(idx,[2 3]) = elem(idx,[3 2]);
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(idx,[2 3]) = bdFlag(idx,[3 2]);
end