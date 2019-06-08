function [node,elem,bdFlag,HB,tree] = bisect(node,elem,markedElem,bdFlag)
%% BISECT bisect a 2-D triangulation.
% 
% [node,elem] = BISECT(node,elem,markedElem) refine the current
% triangulation by bisecting marked elements and minimal neighboring
% elements to get a conforming and shape regular triangulation. Newest
% vertex bisection is implemented. markedElem is a vector containing the
% indices of elements to be bisected. It could be a logical vector of
% length size(elem,1). 
% 
% [node,elem,bdFlag] = BISECT(node,elem,markedElem,bdFlag) returns the
% updated bdFlag after the bisection. It will be used for PDEs with mixed
% boundary conditions.
% 
% [node,elem,bdFlag,HB,tree] = BISECT(node,elem,markedElem,bdFlag)
% returns HB and tree arrays.
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
%      [node,elem] = bisect(node,elem,'all');
%      figure(1); subplot(1,3,2); showmesh(node,elem);
%      bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
%      [node,elem,bdFlag] = bisect(node,elem,[1 4],bdFlag);
%      figure(1); subplot(1,3,3); showmesh(node,elem);
%
% See also bisect3, coarsen, coarsen3, nodeinterpolate, eleminterpolate.
%
% Reference page in Help browser
% <a href="matlab:ifem meshdoc">ifem meshdoc</a>
% <a href="matlab:ifem bisectdoc">ifem bisectdoc</a> 

% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Set up
HB = []; tree = []; 
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('markedElem','var'), markedElem = (1:size(elem,1))'; end
if isempty(markedElem), return; end
if strcmp(markedElem,'all'), markedElem = (1:size(elem,1))'; end
if islogical(markedElem), markedElem = find(markedElem); end

%% Construct auxiliary data structure
T = auxstructure(elem);
neighbor = T.neighbor; elem2edge = T.elem2edge; edge = T.edge;
clear T;
%[neighbor,elem2edge,edge] = auxstructurec(int32(elem));
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Add new nodes
isCutEdge = false(NE,1);
while sum(markedElem)>0
    isCutEdge(elem2edge(markedElem,1)) = true;
    refineNeighbor = neighbor(markedElem,1);
    markedElem = refineNeighbor(~isCutEdge(elem2edge(refineNeighbor,1)));
end
edge2newNode = zeros(NE,1,'uint32');
edge2newNode(isCutEdge) = N+1:N+sum(isCutEdge);
HB = zeros(sum(isCutEdge),3,'uint32');
HB(:,1) = edge2newNode(isCutEdge);
HB(:,[2 3]) = edge(isCutEdge,[1 2]);
node(HB(:,1),:) = (node(HB(:,2),:) + node(HB(:,3),:))/2;

%% Refine marked elements
Nb = 0; tree = zeros(3*NT,3,'uint32');
for k = 1:2
    t = find(edge2newNode(elem2edge(:,1))>0);
    newNT = length(t);
    if (newNT == 0), break; end
    L = t; R = NT+1:NT+newNT;
    p1 = elem(t,1); p2 = elem(t,2); p3 = elem(t,3);
    p4 = edge2newNode(elem2edge(t,1));
    elem(L,:) = [p4, p1, p2];
    elem(R,:) = [p4, p3, p1];
	if nargin==4 && ~isempty(bdFlag) % Refine boundary edges
   		bdFlag(R,[1 3]) = bdFlag(t,[2 1]);
   		bdFlag(L,[1 2]) = bdFlag(t,[3 1]);
        bdFlag(L,3) = 0;
    else
        bdFlag = [];
	end
    tree(Nb+1:Nb+newNT,1) = L;
    tree(Nb+1:Nb+newNT,2) = L;
    tree(Nb+1:Nb+newNT,3) = R;
    elem2edge(L,1) = elem2edge(t,3);
    elem2edge(R,1) = elem2edge(t,2);
    NT = NT + newNT; Nb = Nb + newNT;
end
tree = tree(1:Nb,:);