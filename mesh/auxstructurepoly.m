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
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~iscell(elem)
    warning('elem:incorrectType',...
       'Input being %s may yield incorrect result, please use cell array',class(elem))
    elem = mat2cell(elem,2);
end
%%
elemVertexNumber = cellfun('length',elem);
maxNv = max(elemVertexNumber);
minNv = min(elemVertexNumber);
NT = size(elem,1);
NallEdge = sum(elemVertexNumber); % # of edge = # of vertices locally in 2D
allEdgeCell = cell(NT,1); % cell array of the global indices of all edges
IdxLocalEdge = cell(NT,1); % the local indexing of the edges
IdxGlobalElem = cell(NT,1); % repmat of global indexing of the elements in cells
idxNv = cell(maxNv,1);
for Nv = minNv:maxNv
    idxNv{Nv} = (elemVertexNumber == Nv);
    elemNv = cell2mat(elem(idxNv{Nv}));
    % number of elements sharing the same # of vertices
    NelemNv = sum(idxNv{Nv}); 
    locEdge = [1:Nv; circshift(1:Nv,-1)]';
    allEdgeNv = zeros(NelemNv*Nv,2);
    for i = 1:Nv
        allEdgeNv(i:Nv:(NelemNv-1)*Nv+i,:) = elemNv(:,locEdge(i,:));
    end
    allEdgeNv = permute(reshape(allEdgeNv',[2, Nv, NelemNv]),[2 1 3]);
    allEdgeCell(idxNv{Nv}) = squeeze(num2cell(allEdgeNv, [1,2]));
    locIdxFacetNv = repmat(1:Nv, [NelemNv, 1]);
    IdxLocalEdge(idxNv{Nv}) = num2cell(locIdxFacetNv, 2);
    locIdxElemNv = repmat(find(idxNv{Nv}), [1, Nv]);
    IdxGlobalElem(idxNv{Nv}) = num2cell(locIdxElemNv,2);
end

%%
if isrow(allEdgeCell); allEdgeCell = allEdgeCell'; end
allEdge = cell2mat(allEdgeCell);
allEdge = sort(allEdge,2);
IdxLocalEdge = cat(2, IdxLocalEdge{:})';
IdxGlobalElem = cat(2, IdxGlobalElem{:})';
%%

[edge, i2, j] = myunique(allEdge);
elem2edge = mat2cell(uint32(j), elemVertexNumber);
i1(j(NallEdge:-1:1)) = NallEdge:-1:1;
i1 = i1';
t1 = IdxGlobalElem(i1); t2 = IdxGlobalElem(i2);
k1 = IdxLocalEdge(i1); k2 = IdxLocalEdge(i2); 
ix = (i1 ~= i2); 
edge2elem = uint32([t1,t2,k1,k2]);
neighbor = uint32(accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT maxNv]));
neighbor = cellfun(@(x)x(x~=0), num2cell(neighbor,2),'UniformOutput',false);
       
%%
T.edge2elem = edge2elem;
T.elem2edge = elem2edge;
T.neighbor = neighbor;
T.edge = uint32(edge);
