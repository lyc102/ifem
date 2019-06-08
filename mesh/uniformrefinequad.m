function [node,elem,bdFlag] = uniformrefinequad(node,elem,bdFlag)
%% UNIFORMREFINEQUAD uniformly refine a 2-D quad mesh.
%
% [node,elem] = uniformrefinequad(node,elem) divides each quad into 4 small
% quad by connecting the middle points of the opposited edge of every quad.
% 
% Author: Huayi Wei < huayiwei1984@gmail.com>
% Modified by Long Chen. Add update of bdFlag.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
T = auxstructurequad(elem);
edge = double(T.edge);
elem2edge = double(T.elem2edge);
clear T;

N = size(node); NT = size(elem); NE = size(edge,1);

%% Add new nodes
node(N+1:N+NE,:) = (node(edge(:,1),:) + node(edge(:,2),:))/2; % middle pts of edges
node(N+NE+1:N+NE+NT,:) = (node(elem(:,1),:) + node(elem(:,2),:) + ...
                          node(elem(:,3),:) + node(elem(:,4),:))/4; % center of elements 

%% Refine each quad into four quads as follows
% 4  - 7 -  3
% | t4 | t3 |
% 8 -  9 -  6
% | t1 | t2 |
% 1 -  5 -  2
edge2newNode = (N+1:N+NE)';
elem2newNode = (N+NE+1:N+NE+NT)';
t = 1:NT;
p(t,1:4) = elem(t,1:4);
p(t,5:8) = edge2newNode(elem2edge(t,1:4));
p(t,9) = elem2newNode;
elem(t,:) = [p(t,1), p(t,5), p(t,9), p(t,8)];
elem(NT+1:2*NT,:) = [p(t,5), p(t,2), p(t,6), p(t,9)];
elem(2*NT+1:3*NT,:) = [p(t,9), p(t,6), p(t,3), p(t,7)];
elem(3*NT+1:4*NT,:) = [p(t,8), p(t,9), p(t,7), p(t,4)]; 

%% Update boundary edges
if (nargin>=3) && (~isempty(bdFlag))
    bdFlag(3*NT+1:4*NT,1) = 0;  % extend bdFlag to fine grid
    bdFlag(NT+1:2*NT,[1 2]) = bdFlag(t,[1 2]); 
    bdFlag(2*NT+1:3*NT,[2 3]) = bdFlag(t,[2 3]); 
    bdFlag(3*NT+1:4*NT,[3 4]) = bdFlag(t,[3 4]);
    bdFlag(t,[2 3]) = 0;
else
    bdFlag = [];
end