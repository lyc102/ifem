function [elem,diamond,tree,elemIdxMap] = uniformcoarsen(elem)
%% UNIFORMCOARSEN uniform coarsening
%
%   [elem,diamond,tree] = UNIFORMCOARSEN(elem) remove all
%   good-to-coarsen nodes. It is mainly used to get multilevel
%   decomposition in multigrid methods. The input matrix elem stands for
%   the fine mesh and the output one for the coarse mesh, whose nodal
%   indices are shifted and thus different with the fine mesh.
% 
%   Local index of two types of diamonds:
% 
%   interori nodes    4        boundary nodes       4
%                   / | \                         / | \
%                  2--1--3                       2--1--3
%                   \ | /
%                     5
% 
%   Note that diamond(2:5) are nodal indices in the coarse mesh.
% 
%   See also: coarsen, bisect, uniformcoarsen3, mg
% 
% Documentation in Help browser 
% <a href="matlab:ifem coarsendoc">ifem coarsendoc</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Find good-to-coarsen nodes
N = max(elem(:));
NT = size(elem,1);
valence = accumarray(elem(:),ones(3*NT,1),[N 1]);
valenceNew = accumarray(elem(:,1),ones(NT,1), [N 1]); % for newest vertex only
allNodes = (1:N)';
intGoodNode = allNodes((valence==valenceNew) & (valence==4));
bdGoodNode = allNodes((valence==valenceNew) & (valence==2));
NTdead = 2*length(intGoodNode) + length(bdGoodNode);
if (NTdead==0)
    diamond = []; tree = []; elemIdxMap = [];
    return;
else
    diamond = zeros(N,5,'int32');
end

%% Remove interiori good-to-coarsen nodes
t2v = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), 1, NT, N);
% Find stars for good-to-coarsen nodes
[ii,jj] = find(t2v(:,intGoodNode));
if length(jj)<0
    error('good nodes are not correct')
end
nodeStar = reshape(ii,4,length(intGoodNode));
ix = (elem(nodeStar(1,:),3) == elem(nodeStar(4,:),2));
iy = (elem(nodeStar(2,:),2) == elem(nodeStar(3,:),3));
% nodeStar(:,ix & iy)  = nodeStar([1 2 3 4],ix & iy);
nodeStar(:,ix & ~iy) = nodeStar([1 3 2 4],ix & ~iy);
nodeStar(:,~ix & iy) = nodeStar([1 4 2 3],~ix & iy);
% nodeStar(:,~ix & ~iy)= nodeStar([1 4 3 2],~ix & ~iy);
t1 = nodeStar(1,:); 
t2 = nodeStar(2,:); 
t3 = nodeStar(3,:);
t4 = nodeStar(4,:);
p1 = elem(t1,1); 
p2 = elem(t1,3); 
p3 = elem(t2,2); 
p4 = elem(t1,2); 
p5 = elem(t3,2);
elem(t1,:) = [p4 p2 p3]; 
elem(t2,1) = 0;
elem(t3,:) = [p5 p3 p2]; 
elem(t4,1) = 0;
diamond(p1,:) = [p1 p2 p3 p4 p5];

%% Remove boundary good-to-coarsen nodes
% Find stars for good-to-coarsen nodes
[ii,jj] = find(t2v(:,bdGoodNode));
if length(jj)<0
    error('good nodes are not correct')
end
nodeStar = reshape(ii,2,length(bdGoodNode));
t5 = nodeStar(1,:); 
t6 = nodeStar(2,:);
p1 = elem(t5,1); 
p2 = elem(t5,3); 
p3 = elem(t6,2); 
p4 = elem(t5,2);
elem(t5,:) = [p4 p2 p3]; 
elem(t6,1) = 0;
diamond(p1,1:5) = [p1 p2 p3 p4 p4];

%% Record tree structure
tree = zeros(NTdead,2,'uint32');
tree(1:NTdead,1) = [t1'; t3'; t5'];
tree(1:NTdead,2) = [t1'; t3'; t5'];
tree(1:NTdead,3) = [t2'; t4'; t6'];

%% Clean node,elem and shift elem and diamond matrices
diamond((diamond(:,1) == 0),:) = [];
isCoarseElem = true(NT,1);
isCoarseElem(elem(:,1) == 0) = false;
elem((elem(:,1) == 0),:) = [];
isCoarseNode = true(N,1);
isCoarseNode(intGoodNode) = false;
isCoarseNode(bdGoodNode) = false;
indexMap = zeros(N,1);
indexMap(isCoarseNode) = 1:sum(isCoarseNode);
elem = indexMap(elem);
diamond(:,2:5) = indexMap(diamond(:,2:5));
% shift tree
elemIdxMap = 1:NT;
elemIdxMap(isCoarseElem) = 1:size(elem,1);
tree(:,1) = elemIdxMap(tree(:,1));