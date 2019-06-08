function II = node2edgematrix1(node,edge,isBdEdge)
%% NODE2EDGEMATRIX1 matrix to transfer nodal element to the linear edge element
%
% II = node2edgematrix1(node,edge,isBdEdge) returns the sparse matrix II
% which an NE by 3N matrix mapping a vector (with length 3N) of linear
% nodal element to second family linear edge element (with length 2*NE). The
% trace is zero at certain boundary edges given by isBdEdge.
% 
% II = node2edgematrix1(node,edge) returns the corresponding matrix without
% boundary conditions.
%
% This function is used in mgMaxwell to construct transfer operators for
% the HX preconditioner.
%
% See also mgMaxwell, mgMaxwell1, mgMaxwell2, node2edgematrix1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if nargin<=2, isBdEdge = []; end
N = size(node,1);   NE = size(edge,1);
edgeVec = node(edge(:,2),:)-node(edge(:,1),:);
% first family linear element
i = repmat((1:NE)',6,1);
j = double([edge(:); edge(:)+N; edge(:)+2*N]);
s = 0.5*[edgeVec(:,1); edgeVec(:,1); edgeVec(:,2); edgeVec(:,2); ...
         edgeVec(:,3); edgeVec(:,3)];
bdEdge = edge(isBdEdge,:);
isBdNode = false(3*N,1);
isBdNode([bdEdge(:); bdEdge(:)+N; bdEdge(:)+2*N]) = true;
idx = ~(isBdEdge(i) | isBdNode(j)); 
II = sparse(i(idx),j(idx),s(idx),2*NE,3*N);
% second family linear element
i = repmat((NE+1:2*NE)',6,1);
j = double([edge(:); edge(:)+N; edge(:)+2*N]);
s = 0.5*[edgeVec(:,1); -edgeVec(:,1); edgeVec(:,2); -edgeVec(:,2); ...
         edgeVec(:,3); -edgeVec(:,3)];
II = II + sparse(i(idx),j(idx),s(idx),2*NE,3*N);