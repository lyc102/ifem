function grad = gradmatrix(edge,isBdEdge)
%% GRADMATRIX matrix for the gradient of a nodal linear element
%
% grad = gradmatrix(edge) returns the sparse matrix grad which an NE by 3N
% matrix mapping linear nodal element (a vector with length N) to lowest
% order linear edge element (a vector with length NE). The trace is zero at
% certain boundary edges given by isBdEdge.
% 
% grad = gradmatrix(edge) returns the corresponding matrix without boundary
% conditions.
%
% This function is used in mgMaxwell to construct transfer operators for
% the HX preconditioner.
%
% See also mgMaxwell, mgMaxwell1, mgMaxwell2, node2edgematrix1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin<=1, isBdEdge = []; end
NE = size(edge,1); N  = double(max(edge(:)));
i = repmat((1:NE)',2,1);
j = double(edge(:));
s = [-ones(NE,1),ones(NE,1)];
bdEdge = edge(isBdEdge,:);
isBdNode = false(N,1);
isBdNode(bdEdge(:)) = true;
idx = ~(isBdEdge(i) | isBdNode(j)); 
% idx = 1:length(i);
grad = sparse(i(idx),j(idx),s(idx),NE,N);