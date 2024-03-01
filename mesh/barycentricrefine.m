function [node,elem,bdFlag] = barycentricrefine(node,elem,bdFlag)
% BARYCENTRICREFINE  barycentric refine a 2-D triangulatio.
%
% [node,elem] = barycentricrefine(node,elem) divides each triangle into 3
% small triangles by adding the barycenter.
%
% Example
%
%      node = [0,0; 1,0; 1,1; 0,1];
%      elem = [2,3,1; 4,1,3];
%      figure(1); subplot(1,3,1); showmesh(node,elem);
%      [node,elem] = uniformrefine(node,elem);
%      figure(1); subplot(1,3,2); showmesh(node,elem);
%      bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
%      [bdNode,bdEdge,isBdNode] = findboundary(elem,bdFlag); 
%      findnode(node,bdNode);
%      findedge(node,bdEdge,'all','draw');
%      [node,elem,bdFlag] = barycentricrefine(node,elem,bdFlag);
%      figure(1); subplot(1,3,3); showmesh(node,elem); hold on
%      [bdNode,bdEdge,isBdNode] = findboundary(elem,bdFlag); 
%      findnode(node,bdNode);
%      findedge(node,bdEdge,'all','draw');
%
% See also uniformrefine, uniformbisect, bisect, barycentricrefine3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Add barycenters to node
N = size(node,1); NT = size(elem,1); 
node(N+1:N+NT,:) = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;

%% Add triangles
%        3
%      / | \
%     /  4  \
%    /  /  \ \
%   1 ------- 2 
t= 1:NT; t1 = t; t2 = NT+t; t3 = 2*NT+t;
p(t,1:3) = elem(t,1:3);
p(t,4) = N+t;
elem(t3,:) = [p(t,1), p(t,2), p(t,4)];
elem(t1,:) = [p(t,2), p(t,3), p(t,4)];
elem(t2,:) = [p(t,3), p(t,1), p(t,4)];

%% Update bdFlag
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlagC = bdFlag;
    bdFlag = zeros(3*NT,3,'uint8');
    bdFlag(t1,3) = bdFlagC(t,1);
    bdFlag(t2,3) = bdFlagC(t,2);
    bdFlag(t3,3) = bdFlagC(t,3);
else
    bdFlag = [];
end