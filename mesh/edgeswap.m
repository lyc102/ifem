function [elem,flag] = edgeswap(node,elem)
%% EDGESWAP swap edges
%
% elem = edgeswap(node,elem) swap diagonals of a convex quadlateral
% formed by two triangles of the mesh (node,elem) such that the summation
% of opposite angles are less than or equal to 180 degree, i.e., it is locally
% Delaunay.
%
% Since the connectivity may be changed, some elementwise quantity, e.g.,
% the material property, if any, should be updated accordingly.
%
% Example:
%  load airfoilperturbmesh
%  meshquality(node,elem);
%  elem = edgeswap(node,elem);
%  meshquality(node,elem);
%
% See also  bdsmoothing, meshsmoothing, rmisopoint
% 
% Author: Huayi Wei <huayiwei1984@gmail.com>
%         Long Chen. help comments and vectorize the part: find independent edges.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


NT = size(elem,1);

%% Compute the three angles of all triagnles
% edge vector
v32 = node(elem(:,3),:) - node(elem(:,2),:);
v13 = node(elem(:,1),:) - node(elem(:,3),:);
v21 = node(elem(:,2),:) - node(elem(:,1),:);
% edge length
s1 = sqrt(sum(v32.^2,2));
s2 = sqrt(sum(v13.^2,2));
s3 = sqrt(sum(v21.^2,2));
% angles
a1 = acosd((s2.^2 + s3.^2 - s1.^2)./(2*s2.*s3));
a2 = acosd((s3.^2 + s1.^2 - s2.^2)./(2*s3.*s1));
a3 = acosd((s1.^2 + s2.^2 - s3.^2)./(2*s1.*s2));
a = [a1,a2,a3];

%% Construct necessary data structure
T = auxstructure(elem);
edge2elem = double(T.edge2elem);
elem2edge = double(T.elem2edge);
clear T;

%% Find non-Delaunay edges
idx1 = edge2elem(:,1) + NT*(edge2elem(:,3)-1);
idx2 = edge2elem(:,2) + NT*(edge2elem(:,4)-1);
asum = a(idx1) + a(idx2);
% Find edges with the summation of opposite angles are bigger than 180
% degree and the edge is not a boundary edge
isNeedEdgeSwap = (asum > 180) & (edge2elem(:,1)~= edge2elem(:,2)); 
if ~any(isNeedEdgeSwap)
   flag = 0;
   return 
end

%% Find independent set of swap edges
% Check that if there exsit some elements which have more than one edge to
% swap. Choose the edge with the biggest summation angle to swap.
isCheckElem = (isNeedEdgeSwap(elem2edge(:,1)) + isNeedEdgeSwap(elem2edge(:,2)) ...
          + isNeedEdgeSwap(elem2edge(:,3)) > 1);
if any(isCheckElem)      
    acheck = asum(elem2edge(isCheckElem,:));
    isNeedEdgeSwap(elem2edge(isCheckElem,:)) = false;
    [tempvar,I] = max(acheck,[],2); %#ok<*ASGLU>
    isNeedEdgeSwap(elem2edge(isCheckElem,I)) = true;
end

%% Swap edges
if any(isNeedEdgeSwap)
    %   3 - 2        3 - 2
    %  / / /   to   / \ / 
    % 4 - 1        4 - 1
    localMap = [2;3;1];
    p1 = edge2elem(isNeedEdgeSwap,1) + NT*(edge2elem(isNeedEdgeSwap,3) - 1);
    p2 = edge2elem(isNeedEdgeSwap,1) + NT*(localMap(edge2elem(isNeedEdgeSwap,3)) - 1);
    p3 = edge2elem(isNeedEdgeSwap,2) + NT*(edge2elem(isNeedEdgeSwap,4) - 1);
    p4 = edge2elem(isNeedEdgeSwap,2) + NT*(localMap(edge2elem(isNeedEdgeSwap,4)) - 1);

    leftElem = [ elem(p1), elem(p3), elem(p4)];
    rightElem = [ elem(p3), elem(p1), elem(p2)];
    elem(edge2elem(isNeedEdgeSwap,1),:) = leftElem;
    elem(edge2elem(isNeedEdgeSwap,2),:) = rightElem;
    flag = 1;
end