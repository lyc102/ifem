function [node,elem,bdFlag,HB,tree] = uniformbisect(node,elem,bdFlag)
%% UNIFORMBISECT uniformly bisect a 2-D triangulation. 
%
% [node,elem] = UNIFORMBISECT(node,elem) divides each triangle into 4 small
% triangles using newest vertex bisection. Note that dividing only into 2
% chirldren maynot result in a conforming mesh.
%
% [node,elem,bdFlag] = UNIFORMBISECT(node,elem,bdFlag) update boundary
% conditions represented by bdFlag.
%
% [node,elem,bdFlag,HB,tree] = UNIFORMBISECT(node,elem) returns HB and tree
% arrays.
%
% - HB(:,1:3) is a hierarchical basis structure for nodes, where
%   HB(:,1) is the global index of new added nodes, and HB(:,2:3) the 
%   global indices of two parent nodes of new added nodes. HB is usful
%   for the interpolation between two grids; see also nodeinterpolate.
% 
% - tree(:,1:3) stores the binary tree of the coarsening. tree(:,1) is the
%   index of parent element in coarsened mesh and tree(:,2:3) are two
%   children indices in original mesh. tree is useful for the interpolation
%   of elementwise function; see also eleminterpolate.
%
% Example
%
%      node = [0,0; 1,0; 1,1; 0,1];
%      elem = [2,3,1; 4,1,3];
%      figure(1); subplot(1,3,1); showmesh(node,elem);
%      [node,elem] = uniformbisect(node,elem);
%      figure(1); subplot(1,3,2); showmesh(node,elem);
%      bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
%      [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
%      figure(1); subplot(1,3,3); showmesh(node,elem);
%
% See also uniformbisect3, uniformrefine, bisect, bisect3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('bdFlag','var'), bdFlag = []; end

%% Construct data structure
totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge, i2, j] = myunique(totalEdge); %#ok<*ASGLU>
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
elem2edge = uint32(reshape(j,NT,3));

%% Add new nodes: middle points of all edges
node(N+1:N+NE,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2;
HB(:,[1 2 3]) = [(N+1:N+NE)', edge(:,1:2)]; 
edge2newNode = uint32((N+1:N+NE)');

%% Bisect each triangle into four triangles
%      1                      1
%    / | \                  / | \
%   6  |  5                /  |  \
%  / \ | / \              /   |   \
% 2 -- 4 -- 3            2 -- 4 -- 3
Nb = 0; tree = zeros(3*NT,2,'uint32');
for k = 1:2 % bisect twice
    t = 1:NT;
    p1 = elem(t,1);
    p2 = elem(t,2);
    p3 = elem(t,3);
    p4 = edge2newNode(elem2edge(t,1));
    elem(t,:) = [p4, p1, p2];
    elem(NT+1:NT+NT,:) = [p4, p3, p1];
    if ~isempty(bdFlag)
    	bdFlag(NT+1:NT+NT,[1 3]) = bdFlag(t,[2 1]);
   		bdFlag(t,[1 2]) = bdFlag(t,[3 1]); 
        bdFlag(t,3) = 0;
    else
        bdFlag = [];
    end   
    tree(Nb+1:Nb+NT,1) = t;
    tree(Nb+1:Nb+NT,2) = t;
    tree(Nb+1:Nb+NT,3) = NT+t;
    elem2edge(t,1) = elem2edge(t,3);
    elem2edge(NT+1:NT+NT,1) = elem2edge(t,2);
    Nb = Nb + NT;
    NT = NT + NT;
end