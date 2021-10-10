function [curlu,volume,curlPhi] = curlu3(node,elem,varargin)
%% CURLU3 curl of a lowest order Nedelec or linear Nedelec function in 3-D.
%
% curlu = curlu3(node,elem,u) compute the gradient of an ND0 or a linear ND
% vector field  u on a 3-D mesh representing by edge basis 
% lambda_i Dlambda_j - lambda_j Dlambda_i, lambda_i Dlambda_j + lambda_j Dlambda_i.
%
% curlu(t, i, e): the i-th component (i=1:3) of the e-th edge (e=1:6) in
% elem(t,:).
% 
% [curlu,volume,curlPhi] = curlu3(node,elem,u) also outputs volume and curlPhi
% which is the curl of ND0 basis (additional part of linear ND is curl free).
%
% [curlPhi,volume] = curlu3(node,elem) if u is not given, the function only
% outputs the volume and curlPhi the curl of ND0 basis.
% 
% Local node index in a tetrahedron
%    _1_
%   / | \
%  2_ | _3
%    -4-
% 6 edges are indexed with incremental fashion.
% The global indexing of the additional curl free bases is assumed to be (NE+1:2*NE).
% 
% See also gradu3, gradbasis3, Maxwell1saddle
%
% Added by Shuhao Cao, Apr 2020. 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Sort elem to ascend ordering
elem = sort(elem,2);
[elem2edge, edge] = dof3edge(elem);
NT = size(elem,1);
[volume,elemSign] = simplexvolume(node,elem);
NE = size(edge,1);
%% Edge vector
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
elem2ve = zeros(NT,3,6);
for e = 1:6 % six edges
   elem2ve(:,:,e) = node(elem(:,locEdge(e,2)),:)-node(elem(:,locEdge(e,1)),:);
end

%% curlPhi is the scaled edge vector
curlPhi(:,:,6) = elem2ve(:,:,1);
curlPhi(:,:,1) = elem2ve(:,:,6);
curlPhi(:,:,2) = -elem2ve(:,:,5);
curlPhi(:,:,3) = elem2ve(:,:,4);
curlPhi(:,:,4) = elem2ve(:,:,3);
curlPhi(:,:,5) = -elem2ve(:,:,2);
curlPhi = curlPhi./repmat(volume.*elemSign,[1 3 6])/3;
% the above ordering is given for positive orientation and should be
% corrected by multiplying the elemSign

%% curl u is the composition of u and curlPhi
if nargin > 2
    u = varargin{1}; 
    if length(u) == 2*NE; u = u(1:NE); end
    curlu = dot(permute(repmat(u(elem2edge),[1,1,3]),[1,3,2]),curlPhi,3);
else
    curlu = curlPhi;
end