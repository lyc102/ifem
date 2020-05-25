function T = auxstructurepoly(elem)
%% AUXSTRUCTUREPOLY auxiliary structure for a 2-D polygonal mesh
%
%  T = auxstructurepoly(elem) constucts the indices map between elements, edges 
%  and nodes, and the boundary information. T is a structure. elem is a
%  cell array, the implementation keeps most Long's original naming
%  traditions unchanged while using cellfun for most parts.
%
%  T.neighbor{1:NT}: cell array of the indices map of neighbor information of elements, 
%  where neighbor(t,i) is the global index of the element on the other side of the 
%  i-th (local) edge of the t-th (global) element. 
%
%  T.elem2edge{1:NT}: cell array of the indices map from elements to edges, 
%  elem2edge{t}(i) is the the global index of the i-th edge of t-th element.
%
%  T.edge(1:NE,1:2): all edges, where edge(e,i) is the global index of the 
%  i-th vertex of the e-th edge, and edge(e,1) < edge(e,2) 
%
%  T.edge2elem(1:NE,1:4): the indices map from edge to element, where 
%  edge2elem(e,1:2) are the global indexes of two elements sharing the e-th
%  edge, and edge2elem(e,3:4) are the local indices of e to edge2elem(e,1:2).
%
%  T.bdEdge(1:Nbd,1:2): boundary edges with counterclockwise oritentation, where
%  bdEdge(e,i) is the global index of the i-th vertex of the e-th edge for
%  i=1,2. The counterclockwise oritentation means that the interior of the domain
%  is on the right moving from bdEdge(e,1) to bdEdge(e,2). 
%
%  To save space all the data type in T is uint32. When use them as a input
%  of sparse(i,j,s,m,n), please change them into double type.
%
%  See auxstructure, auxstructure3, auxstructurequad.
%
%  Notes on the vectorization used: https://scaomath.github.io/blog/matlab-cell-vectorization/
%
%  Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~iscell(elem)
    warning('elem:incorrectType',...
       'Input being %s may yield incorrect result, please use cell array',class(elem))
    elem = num2cell(elem,2);
end
%%
elemVertNum = cellfun('length',elem);
maxNv = max(elemVertNum);
minNv = min(elemVertNum);
NT = size(elem,1);
NallEdge = sum(elemVertNum); % # of edge = # of vertices locally in 2D
allEdgeCell = cell(NT,1); % cell array of the global indices of all edges
IdxLocalEdge = cell(NT,1); % the local indexing of the edges
IdxGlobalElem = cell(NT,1); % repmat of global indexing of the elements in cells
for nV = minNv:maxNv
    isNv = (elemVertNum == nV);
    if ~any(isNv); continue; end
    
    elemNv = cell2mat(elem(isNv));
    % number of elements sharing the same # of vertices
    NelemNv = sum(isNv);
    
    %% new implementation modified from PoissonVEM routine
    elemNvShift = circshift(elemNv,[0,-1])';
    elemNv = elemNv';
    allEdgeCell(isNv) = mat2cell([elemNv(:), elemNvShift(:)], nV*ones(1,NelemNv));
    
    %% recording local indexing for edge2elem
    locIdxFacetNv = repmat(1:nV, [NelemNv, 1]);
    IdxLocalEdge(isNv) = num2cell(locIdxFacetNv, 2);
    locIdxElemNv = repmat(find(isNv), [1, nV]);
    IdxGlobalElem(isNv) = num2cell(locIdxElemNv,2);
end

%%
if isrow(allEdgeCell); allEdgeCell = allEdgeCell'; end
allEdge = cell2mat(allEdgeCell);
IdxLocalEdge = [IdxLocalEdge{:}]';
IdxGlobalElem = [IdxGlobalElem{:}]';
%%
[edge, i2, j] = myunique(sort(allEdge,2));
elem2edge = mat2cell(uint32(j), elemVertNum);
i1(j(NallEdge:-1:1)) = NallEdge:-1:1;
i1 = i1';
t1 = IdxGlobalElem(i1); t2 = IdxGlobalElem(i2);
k1 = IdxLocalEdge(i1); k2 = IdxLocalEdge(i2); 
ix = (i1 ~= i2); 
edge2elem = uint32([t1,t2,k1,k2]);
neighbor = uint32(accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT maxNv]));
neighbor = cellfun(@(x)x(x~=0), num2cell(neighbor,2),'UniformOutput',false);
bdEdge = allEdge(i1(i1 == i2), :);

%%
T.edge2elem = edge2elem;
T.elem2edge = elem2edge;
T.neighbor = neighbor;
T.edge = uint32(edge);
T.allEdge = uint32(allEdge);
T.bdEdge = uint32(bdEdge);
