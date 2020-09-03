function [elem,idx,volume,bdFlag] = fixorder3(node,elem,bdFlag)
%% FIXORDER3 fix orientation of tetrahedron 
% 
%   elem = FIXORDER3(node,elem) computes signed volume of all tetrahedron
%   in the triangulation and switch the vertices such that all signed
%   volume is positive.
%   
%   [elem,idx,volume] = FIXORDER3(node,elem) also outputs the index set of
%   elements whose volume is negative.
%
%   [elem,idx,volume,bdFlag] = FIXORDER3(node,elem,bdFlag) changes the bdFlag
%   for boundary conditions. 
%
% See also fixorder
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin==2)
    bdFlag = [];
end
% compute volume and elemSign of each tetrahedron
[volume,elemSign] = simplexvolume(node,elem);
% find tetrahedron with negative volume and switch the vertices
idx = find(elemSign==-1); 
elem(idx,[2 3]) = elem(idx,[3 2]);
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(idx,[2 3]) = bdFlag(idx,[3 2]);
end