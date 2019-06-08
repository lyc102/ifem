function [elem,idx,volume] = fixorder3(node,elem)
%% FIXORDER3 fix orientation of tetrahedron 
% 
%   elem = FIXORDER3(node,elem) computes signed volume of all tetrahedron
%   in the triangulation and switch the vertices such that all signed
%   volume is positive.
%   
%   [elem,idx,volume] = FIXORDER3(node,elem) also outputs the index set of
%   elements whose volume is negative.
%
% See also fixorder
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

% compute volume and elemSign of each tetrahedron
[volume,elemSign] = simplexvolume(node,elem);
% find tetrahedron with negative volume and switch the vertices
idx = find(elemSign==-1); 
elem(idx,[2 3]) = elem(idx,[3 2]);