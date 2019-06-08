function [node,elem] = perturbmesh(node,elem,theta)
% PERTURBMESH perturb all the interior nodes with the local size.
%   
%   [node,elem] = perturbmesh(node,elem,theta), where node and elem is the
%   mesh data, theta is the perturbution parameter.
%
%  See also
%
%  TODO: .
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(node,1);
bdNodeIdx = findboundary(elem);
area = abs(simplexvolume(node,elem));
patchArea = sparse(elem,ones(size(elem,1),3),area*[1,1,1],N,1);
h = sqrt(patchArea/6);
dp = theta*[h h].*(2*rand(N,2)-1);
oldNode = node(bdNodeIdx,:);
node = node + dp;
node(bdNodeIdx,:) = oldNode;