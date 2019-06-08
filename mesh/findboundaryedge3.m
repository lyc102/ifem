function [bdEdge,bdNode,isBdEdge,isBdNode] = findboundaryedge3(edge,elem2edge,bdFlag)
%% FINDBOUNDARYEDGE3 find boundary edges of a three dimensional mesh
%
% [bdEdge,bdNode] = findboundaryedge3(edge,elem2edge,bdFlag) finds all boundary
% edges and boundary nodes of 3D mesh using bdFlag.
% 
% For a 2D mesh, use findboundary
%
% See also findboundary3, findboundary
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

% Find boundary edges and nodes
NE = size(edge,1);
N = max(edge(:));
isBdEdge = false(NE,1);
isBdNode = false(N,1);
isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
bdEdge = edge(isBdEdge,:);
isBdNode(bdEdge) = true;
bdNode = find(isBdNode);