function [node,elem,bdFlag,HB,indexMap,tree]= ...
         coarsen3(node,elem,markedElem,bdFlag,HB)
%% COARSEN3 coarse a 3-D triangulation.
%
% [node,elem,HB] = COARSEN3(node,elem,markedElem,HB) removes good-to-coraen
% nodes whose star are marked for coarsening. Unlike the 2-D version
% carsen, one additional data structure HB is needed. The local index of HB
% is 2 -- 1 -- 3, i.e., HB(:,1) is the middle point of the edge formed by
% HB(:,2:3). HB(:,4) is used to store the generation of the node HB(:,1).
% Note that it is assumed elem(:,4) records the newest vertex added during
% the bisection; see bisect3.
% 
% [node,elem,HB,bdFlag] = COARSEN3(node,elem,markedElem,HB,bdFlag) updates
% boundary conditions represented by bdFlag.
%
% [node,elem,HB,bdFlag,indexMap,tree] = COARSEN3(node,elem,markedElem,HB,bdFlag)
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
%
% Example
%
%     [node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2);
%     bdFlag = setboundary3(node,elem,'Dirichlet');
%     [node,elem,bdFlag,HB] = bisect3(node,elem,1,bdFlag,HB);
%     figure(2); subplot(1,2,1);
%     showmesh3(node,elem,[],'FaceAlpha',0.15); view([210 8]);
%     findnode3(node);
%     [node,elem,bdFlag,HB,indexMap,tree] = coarsen3(node,elem,'all',bdFlag,HB);
%     figure(2); subplot(1,2,2);
%     showmesh3(node,elem,[],'FaceAlpha',0.15); view([210 8]);
%
% 
% <a href="matlab:ifem coarsendoc">coarsen3doc</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

tree = []; 
indexMap  = (1:size(node,1))'; % default output
if isempty(markedElem), return; end
if strcmp(markedElem,'all'), markedElem = (1:size(elem,1))'; end

%% Find good-to-coarsen nodes
N = size(node,1); 
NT = size(elem,1);
generation = zeros(N,1);
generation(HB(:,1)) = HB(:,4);
valence = accumarray(elem(:),ones(4*NT,1),[N 1]);
valenceNew = accumarray(elem(:,4),ones(NT,1), [N 1]);
valenceMarked = accumarray(elem(markedElem,4),ones(length(markedElem),1),[N 1]);
isGoodNode = (valence == valenceNew) ... coarsen new added vertex
           & (valenceNew == valenceMarked) ... all elements of new vertex are marked
           & (generation>0); % not coarsen the initial grid.
if ~any(isGoodNode)
    return
end
t = find(isGoodNode(elem(:,4))); % all elements contaning good nodes

%% Coarsen boundary faces
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdElem = t(bdFlag(t,4)~=0);
    if ~isempty(bdElem)
        isBdGoodNode = false(N,1);
        isBdGoodNode(elem(bdElem,4)) = true;
        bdElem = t(isBdGoodNode(elem(t,4)));
        T = auxstructure3(elem(bdElem,:));  % only the structure of bdElem
        leftNode = HB(elem(bdElem,4),2);
        for k = 1:3
            kidx = (elem(bdElem,k)==leftNode); % idx of bdElem
            tl = bdElem(kidx);                 % left element
            tr = bdElem(T.neighbor(kidx,k));   % matched neighbor 
            bdFlag(tl,k) = bdFlag(tr,4);    
        end
    end
else
	bdFlag=[];
end

%% Coarsen elements and record tree
leftNode = HB(elem(t,4),2);
idx = (elem(t,1)==leftNode) | (elem(t,2)==leftNode) | (elem(t,3)==leftNode);
tl = t(idx); 
tr = t(~idx);
if length(tl) ~= length(tr)
    error('The mesh is not obtained by bisect3');
end
Nl = length(tl); % number of refined elements
% sort tl and tr by the common vertices such that tl and tr matches
leftNode = HB(elem(tl,4),2);
rightNode = HB(elem(tr,4),3);
temptl = zeros(Nl,2);   
idx = (elem(tl,1) == leftNode);
temptl(idx,:) = elem(tl(idx),[2 3]);
idx = (elem(tl,2) == leftNode);
temptl(idx,:) = elem(tl(idx),[1 3]);
idx = (elem(tl,3) == leftNode);
temptl(idx,:) = elem(tl(idx),[1 2]);
temptl = sort(temptl,2);
temptr = zeros(Nl,2);
idx = (elem(tr,1) == rightNode);
temptr(idx,:) = elem(tr(idx),[2 3]);
idx = (elem(tr,2) == rightNode);
temptr(idx,:) = elem(tr(idx),[1 3]);
idx = (elem(tr,3) == rightNode);
temptr(idx,:) = elem(tr(idx),[1 2]);
temptr = sort(temptr,2);
[temptl, Il] = sortrows([temptl elem(tl,4)]); %#ok<*ASGLU>
[temptr, Ir] = sortrows([temptr elem(tr,4)]);
tl = tl(Il);
tr = tr(Ir);
% record tree structure
tree = zeros(Nl,3);
tree(:,1) = tl;
tree(:,2) = tl;
tree(:,3) = tr;
% coarsen tl to t and will remove tr
elem(tl,4) = HB(elem(tl,4),3);

%% Sort element nodes by generations
% the newest node is the node with maxmum generation
if (length(tl)==1)
    idx = max(transpose(generation(elem(t1,:)))); 
else
    [tempvar,idx] = max(generation(elem(tl,:)),[],2); 
end
elem(tl((idx==1)),1:4) = elem(tl((idx==1)),[2 4 3 1]);
elem(tl((idx==2)),1:4) = elem(tl((idx==2)),[3 4 1 2]);
elem(tl((idx==3)),1:4) = elem(tl((idx==3)),[4 2 1 3]);
if exist('bdFlag','var') && ~isempty(bdFlag)
    bdFlag(tl((idx==1)),1:4) = bdFlag(tl((idx==1)),[2 4 3 1]);
    bdFlag(tl((idx==2)),1:4) = bdFlag(tl((idx==2)),[3 4 1 2]);
    bdFlag(tl((idx==3)),1:4) = bdFlag(tl((idx==3)),[4 2 1 3]);
end

%% Clean and shift index
elem(tr,:) = [];
inCoarse = true(NT,1);
inCoarse(tr) = false;
elemidxMap = zeros(NT,1);
elemidxMap(inCoarse) = 1:size(elem,1); 
tree(:,1) = elemidxMap(tree(:,1));
if (nargin==5) && ~isempty(bdFlag)
    bdFlag(tr,:) = [];
end
node(isGoodNode,:) = [];
HB(isGoodNode,:) = [];
indexMap = zeros(N,1);
indexMap(~isGoodNode) = 1:size(node,1);
elem = indexMap(elem);
HB(:,1:3) = indexMap(HB(:,1:3));