function [elem2edge,edge,elem2edgeSign,edgeSign] = dofedgepoly(elem)
%% DOFEDGEPOLY dof structure for edges on a polygonal mesh.
%
% [elem2edge,edge,elem2edgeSign,edgeSign] = DOFEDGE(elem) constructs data
% structure for finite elements associated to edges including CR
% nonconforming element, Ravairt-Thomas element, and Nedelec element etc.
%
% In the input elem is the connectivity cell array for a 2-D polygonal mesh. In
% the output
%
% - elem2edge: the elementwise pointer from elem to edge indices. In each
% polygon, elem2edge{t}(j) is the global index of the j-th edge consisting of
% [j j+1]-th vertices of the t-th element, assuming the vertex indexing is cyclic. 
%
% - edge: the edge matrix satisfying edge(:,1)<edge(:,2). The orientation
% of edge is induced by the ordering of vertices, i.e., from the vertex
% with a smaller global index to a bigger one.
%
% - elem2edgeSign: records the consistency of the local and global edge
% orientation. The orientation of local edges is the induced oritentation
% of polygon, i.e., quadrilateral: [1 2; 2 3; 3 4; 4 1]. 
%
% - edgeSign: edgeSign equals 1 if the edge is consistent with the local edge of its
% first neighbor element in edge2elem(:,1), equals -1 otherwise.
%
% See also dofedge, dof3edge
%
% Remark on the linear indexing used for cell array: 
% https://scaomath.github.io/blog/matlab-cell-vectorization/
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%%
if ~iscell(elem)
    warning('elem:incorrectType',...
       'Input being %s may yield incorrect result, please use cell array',class(elem))
    elem = num2cell(elem,2);
end

%% Generate edge and elem2edge
T = auxstructurepoly(elem);
elem2edge = T.elem2edge;
edge = T.edge;
edge2elem = T.edge2elem;
totalEdge = T.allEdge;
NE = size(edge,1);
elemVertNum = cellfun('length',elem);

%% The sign of edges
elem2edgeSign = ones(sum(elemVertNum),1);
idx = (totalEdge(:,1)>totalEdge(:,2));
elem2edgeSign(idx) = -1;
elem2edgeSign = mat2cell(elem2edgeSign,elemVertNum);

%% Consistency of oritentation of edges
% edgeSign equals 1 if the edge is consistent with the local edge of its
% first neighboring element, equals -1 otherwise.
if nargout >3
    edgeSign = -ones(NE,1); 
    elemSub2ind = [elem{edge2elem(:,1)}]';
    % linear indexing for a cell array
    idxEdge = [0; cumsum(elemVertNum(edge2elem(1:end-1,1)))] + double(edge2elem(:,3));
    isConsistentEdge = (elemSub2ind(idxEdge) == edge(:,1)); 
    edgeSign(isConsistentEdge) = 1;
end
