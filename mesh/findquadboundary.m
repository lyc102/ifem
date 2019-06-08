function [bdNode,bdEdge,isBdNode,isBdElem] = findquadboundary(elem)
%% FINDQUADBOUNDARY finds the boundary of a quad mesh
%
%   [bdNode,bdEdge,isBdNode] = findquadboundary(elem) finds boundary nodes and
%   edges for a 2-dimensional quad mesh. Note only the topological structure of
%   the mesh is needed.
%
%   See also 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[NT, NV] = size(elem);
totalEdge = zeros(NV*NT,2);
totalEdge(:,1) = elem(:);
totalEdge(:,2) = elem([NT+1:NV*NT,1:NT]');
totalEdge(:) = sort(totalEdge,2);
[i,j,s] = find(sparse(double(totalEdge(:,1)),double(totalEdge(:,2)),1));
bdEdge = [i(s==1),j(s==1)];
isBdNode = false(max(elem(:)),1); 
isBdNode(bdEdge) = true;
bdNode = find(isBdNode);
isBdElem = isBdNode(elem(:,1)) | isBdNode(elem(:,2)) | isBdNode(elem(:,3)) | isBdNode(elem(:,4));