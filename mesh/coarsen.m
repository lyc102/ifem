function [node,elem,bdFlag,indexMap,tree] = coarsen(node,elem,markedElem,bdFlag)
%% COARSEN coarsen a 2-D triangulation.
%
% [node,elem] = COARSEN(node,elem,markedElem) removes good-to-coarsen
% nodes whose star are marked for coarsening
%
% [node,elem,bdFlag] = COARSEN(node,elem,markedElem,bdFlag) updates
% boundary conditions represented by bdFlag.
%
% [node,elem,bdFlag,indexMap,tree] = COARSEN(node,elem,markedElem,bdFlag)
% outputs two additional information: indexMap and tree. 
%
% - indexMap is the map between nodes in the fine mesh (node in the input)
%   to that in the coarse mesh (node in the output). For example,
%   indexMap(10) = 6 means the 10-th node in the fine grid  is now the 6-th
%   node in the coarse one. indexMap is useful for the interpolation of
%   nodewise function; see also nodeinterpolate
%
% - tree(:,1:3) stores the binary tree of the coarsening. tree(:,1) is the
%   index of parent element in coarsened mesh and tree(:,2:3) are two
%   children indices in original mesh.
%
% Example
%
%     load data/Lshapemesh;
%     set(gcf,'Units','normal'); set(gcf,'Position',[0.25,0.25,0.7,0.4]);
%     subplot(1,3,1); showmesh(node,elem);
%     [node,elem] = coarsen(node,elem,'all');
%     subplot(1,3,2); showmesh(node,elem);
%     [node,elem] = coarsen(node,elem,'all');
%     subplot(1,3,3); showmesh(node,elem);
% 
% Reference page in Help browser
%  <a href="matlab:ifem coarsendoc">ifem coarsendoc</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = max(elem(:)); NT = size(elem,1);
tree = []; indexMap  = (1:size(node,1))';
if ~exist('node','var') || isempty(node), node = []; end
if ~exist('bdFlag','var') || isempty(bdFlag), bdFlag = []; end
if isempty(markedElem), return; end
if strcmp(markedElem,'all'), markedElem = (1:size(elem,1))'; end
if islogical(markedElem), markedElem = find(markedElem); end

%% Find good-to-coarsen nodes
valence = accumarray(elem(:),ones(3*NT,1),[N 1]);
markedVal = accumarray(elem(markedElem,1),ones(length(markedElem),1),[N 1]);
isIntGoodNode = ((markedVal==valence) & (valence==4));
isBdGoodNode = ((markedVal==valence) & (valence==2));
NTdead = 2*sum(isIntGoodNode) + sum(isBdGoodNode); 
if (NTdead == 0), return; end

%% Remove interiori good-to-coarsen nodes
t2v = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), 1, NT, N);
% Find stars for good-to-coarsen nodes
[ii,jj] = find(t2v(:,isIntGoodNode));
if length(jj)<0
    error('number of good nodes are not correct')
end
nodeStar = reshape(ii,4,sum(isIntGoodNode));
isIntNode = false(size(nodeStar,2),1);
% isIntNode is used to exclude those bd nodes whose val = 4
% case: 1 2 3 4
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(4,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(2,:),3));  
nodeStar(:,idx)  = nodeStar([1 2 3 4],idx);
isIntNode(idx) = true;
% case: 1 2 4 3
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(3,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(2,:),3));  
nodeStar(:,idx)  = nodeStar([1 2 4 3],idx); 
isIntNode(idx) = true;
% case: 1 3 2 4
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(4,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(3,:),3));  
nodeStar(:,idx)  = nodeStar([1 3 2 4],idx); 
isIntNode(idx) = true;
% case: 1 3 4 2
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(2,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(3,:),3));  
nodeStar(:,idx)  = nodeStar([1 3 4 2],idx); 
isIntNode(idx) = true;
% case: 1 4 2 3
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(3,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(4,:),3));  
nodeStar(:,idx)  = nodeStar([1 4 2 3],idx); 
isIntNode(idx) = true;
% case: 1 4 3 2
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(2,:),2)) & ...
      (elem(nodeStar(1,:),2) == elem(nodeStar(4,:),3));  
nodeStar(:,idx)  = nodeStar([1 4 3 2],idx); 
isIntNode(idx) = true;
% merge t1 with t2, and t3 with t4
t1 = nodeStar(1,isIntNode); 
t2 = nodeStar(2,isIntNode); 
t3 = nodeStar(3,isIntNode);
t4 = nodeStar(4,isIntNode);
p2 = elem(t1,3); 
p3 = elem(t2,2); 
p4 = elem(t1,2); 
p5 = elem(t3,2);
elem(t1,:) = [p4 p2 p3]; 
elem(t2,1) = 0;
elem(t3,:) = [p5 p3 p2]; 
elem(t4,1) = 0;
% update isIntGoodNode
intGoodNode = find(isIntGoodNode);
isIntGoodNode(intGoodNode(~isIntNode)) = false;

%% Remove boundary good-to-coarsen nodes
% Find stars for good-to-coarsen nodes
[ii,jj] = find(t2v(:,isBdGoodNode));
if length(jj)<0
    error('number of good nodes are not correct')
end
nodeStar = reshape(ii,2,sum(isBdGoodNode));
idx = (elem(nodeStar(1,:),3) == elem(nodeStar(2,:),2));
nodeStar(:,idx)  = nodeStar([2 1],idx); 
t5 = nodeStar(1,:); 
t6 = nodeStar(2,:);
p1 = elem(t5,1);
p2 = elem(t5,3); 
p3 = elem(t6,2); 
p4 = elem(t5,2);
if ~isempty(node)
    v13 = node(p1,:) - node(p3,:);
    v12 = node(p1,:) - node(p2,:);
    v23 = node(p2,:) - node(p3,:);
    % Corner/feature points could be good nodes. Remove these points will
    % change the shape of the domain. We check if nodes are feature points by
    % computing the length differences.
    lengthDiff = (sum(v12.^2,2) + sum(v13.^2,2) - sum(v23.^2,2))./sum(v23.^2,2);
    idx = (sqrt(lengthDiff) < 1e-3);
else
    idx = true(length(t5),1);
end
elem(t5(idx),:) = [p4 p2 p3]; 
elem(t6(idx),1) = 0;
bdGoodNode = find(isBdGoodNode);
isBdGoodNode(bdGoodNode(~idx)) = false;

%% Update boundary edges
if (nargin==4) && (~isempty(bdFlag))
	bdFlag(t1,:) = [bdFlag(t1,2) bdFlag(t2,1) bdFlag(t1,1)];
	bdFlag(t3,:) = [bdFlag(t3,2) bdFlag(t4,1) bdFlag(t3,1)];	
	bdFlag(t5,:) = [bdFlag(t5,2) bdFlag(t6,1) bdFlag(t5,1)];
    bdFlag((elem(:,1) == 0),:) = [];
else
	bdFlag=[];
end

%% Record tree structure
NTdead = 2*sum(isIntGoodNode) + sum(isBdGoodNode);
tree = zeros(NTdead,3,'uint32');
tree(1:NTdead,1) = [t1'; t3'; t5(idx)'];
tree(1:NTdead,2) = [t1'; t3'; t5(idx)'];
tree(1:NTdead,3) = [t2'; t4'; t6(idx)'];

%% Clean node and elem matrices
isRemoved = (elem(:,1) == 0);
elem(isRemoved,:) = [];
inCoarse = true(NT,1);
inCoarse(isRemoved) = false;
elemidxMap = zeros(NT,1);
elemidxMap(inCoarse) = 1:size(elem,1); 
tree(:,1) = elemidxMap(tree(:,1));
if ~isempty(node)
    node(isIntGoodNode | isBdGoodNode,:) = [];
end
indexMap = zeros(N,1);
Ndead = sum(isIntGoodNode) + sum(isBdGoodNode);
indexMap(~(isIntGoodNode | isBdGoodNode)) = 1:(N-Ndead);
elem = indexMap(elem);