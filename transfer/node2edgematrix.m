function II = node2edgematrix(node,edge,isBdEdge)
%% NODE2EDGEMATRIX matrix to transfer nodal element to the lowest order edge element
%
% II = node2edgematrix(node,edge,isBdEdge) returns the sparse matrix II
% which is an NE by 3N matrix mapping a vector (with length 3N) of linear
% nodal element to lowest order linear edge element (with length NE). The
% trace is zero at certain boundary edges given by isBdEdge.
% 
% II = node2edgematrix(node,edge) returns the corresponding matrix without
% boundary conditions.
%
% This function is used in mgMaxwell to construct transfer operators for
% the HX preconditioner.
%
% See also mgMaxwell, mgMaxwell1, mgMaxwell2, node2edgematrix1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

dim = size(node,2);  % dimension
N = size(node,1);   NE = size(edge,1);
if nargin<=2, isBdEdge = []; end
if isempty(isBdEdge), isBdEdge = false(NE,1); end
edgeVec = node(edge(:,2),:)-node(edge(:,1),:);
if dim == 2
    i = repmat((1:NE)',4,1);
    j = double([edge(:); edge(:) + N]);
    s = 0.5*[edgeVec(:,1); edgeVec(:,1); edgeVec(:,2); edgeVec(:,2)];
elseif dim == 3
    i = repmat((1:NE)',6,1);
    j = double([edge(:); edge(:) + N; edge(:) + 2*N]);
    s = 0.5*[edgeVec(:,1); edgeVec(:,1); edgeVec(:,2); edgeVec(:,2); ...
             edgeVec(:,3); edgeVec(:,3)];
end
bdEdge = edge(isBdEdge,:);
isBdNode = false(dim*N,1);
if dim == 2
    isBdNode([bdEdge(:); bdEdge(:) + N]) = true;
elseif dim == 3
    isBdNode([bdEdge(:); bdEdge(:) + N; bdEdge(:) + 2*N]) = true;
end
idx = ~(isBdEdge(i) | isBdNode(j)); 
% idx = 1:length(i);
II = sparse(i(idx),j(idx),s(idx),NE,dim*N);