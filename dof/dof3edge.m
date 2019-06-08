function [elem2edge,edge,elem2edgeSign] = dof3edge(elem)
%% DOF3EDGE dof structure for the lowest order edge elements.
%
% [elem2edge,edge] = DOF3EDGE(elem) constructs data structure
% for the lowest order edge element. elem is the connectivity matrix for a
% 3-D triangulation. elem2edge is the elementwise pointer from elem to dof
% indices. edge is the edge matrix. The orientation of edge is
% from the node with smaller index to bigger one. The orientation of local
% edges of a tetrahedron formed by [1 2 3 4] is lexicographic order:
%           [1 2], [1 3], [1 4], [2 3], [2 4], [3 4]
%
% [elem2edge,edge,elem2edgeSign] = DOF3EDGE(elem) also output elem2edgeSign
% which records the consistency of the local and global edge orientation.
% If elem is ascend ordered, then elem2edgeSign is 1 and not needed.
%
% See also dofedge.
%
% Doc: <a href="matlab:ifem dof3edgedoc">dof3edgedoc</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

NT = size(elem,1);
totalEdge = int32([elem(:,[1 2]); elem(:,[1 3]); elem(:,[1 4]); ...
                   elem(:,[2 3]); elem(:,[2 4]); elem(:,[3 4])]);
sortedTotalEdge = sort(totalEdge,2);
[edge,i2,j] = myunique(sortedTotalEdge); %#ok<*ASGLU>
elem2edge = uint32(reshape(j,NT,6));
direction = ones(6*NT,1);
idx = (totalEdge(:,1)>totalEdge(:,2));
direction(idx) = -1;
elem2edgeSign = reshape(direction,NT,6);