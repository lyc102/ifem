function [node,elem,bdNode] = bdsmoothing(node,elem,step)
%% BDSMOOTHING improves the boundary of a triangular mesh
%
% [node,elem] = BDSMOOTHING(node,elem) improves the boundary of a
% triangular mesh by moving boundary nodes while keeping the shape of the
% domain unchanged. Therefore corner vertices are fixed. Here a vertex is
% called a corner vertex if it doesn't lie on the line segement connecting
% its two neighboring vertices. A nice feature of this function is that
% users do not need to provide boundary vertices and corner vertices. The
% function will find all boundary vertices and determine which are corners
% by using topological and geometrical information provided by (node,elem).
%
% [node,elem,bdNode] = BDSMOOTHING(node,elem,k) performs k steps smoothing
% and return boundary vertices bdNode. The default setting is k=2. 
%
% Example:
%  load airfoilperturb
%  subplot(2,2,1); showmesh(node,elem);
%  subplot(2,2,2); meshquality(node,elem);
%  [node,elem] = bdsmoothing(node,elem);
%  subplot(2,2,3); showmesh(node,elem);
%  subplot(2,2,4); meshquality(node,elem);
%
% See also meshsmoothing, edgeswap, rmisopoint
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('step','var'), step =2; end
%% Find boundary edges with positive orientation
T = auxstructure(elem);
bdEdge = T.bdEdge; 
clear T;

%% Find boundary nodes and corner nodes
isbdNode(bdEdge(:)) = true;
bdNode = find(isbdNode');
N = size(node,1);   NE = size(bdEdge,1);
v2e = zeros(N,2);
v2e(bdEdge(:,2),1) = 1:NE;
v2e(bdEdge(:,1),2) = 1:NE;
neighbor(bdNode,2) = bdEdge(v2e(bdNode,2),2);
neighbor(bdNode,1) = bdEdge(v2e(bdNode,1),1);
l1 = sqrt(sum((node(bdNode,:) - node(neighbor(bdNode,1),:)).^2,2));
l2 = sqrt(sum((node(neighbor(bdNode,2),:) - node(bdNode,:)).^2,2));
l3 = sqrt(sum((node(neighbor(bdNode,2),:) - node(neighbor(bdNode,1),:)).^2,2));
corner = bdNode((l1+l2-l3)./l3>1e-6);

%% Smoothing boundary nodes
for k = 1:step
    mid = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;
    force1 = mid - node(bdEdge(:,1),:);
    force2 = mid - node(bdEdge(:,2),:);
    dE = zeros(N,2);
    dE(:,1) = accumarray(bdEdge(:),[force1(:,1);force2(:,1)],[N 1]);
    dE(:,2) = accumarray(bdEdge(:),[force1(:,2);force2(:,2)],[N 1]);
    dE(corner,:) = 0;
    node(bdNode,:) = node(bdNode,:) + dE(bdNode,:);
end