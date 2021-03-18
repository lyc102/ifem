function curlu = curlu3CR(elem2face,u,Dlambda)
%% CURLU3CR curl of an C-R vector field in 3-D.
%
% curlu = curlu3CR(elem2face,u,Dlambda) compute the 
% curl of an C-R vector field u on a 3-D mesh representing by 
% face bases in each component. u should be of an NF by 3 array.
% 
% Local node index in a tetrahedron, the face's indexing follows the node
% opposite to the face of interest
%    _1_
%   / | \
%  2_ | _3
%    -4-
% 
% See also curlu3
% 
% Reference: curlu3CR is used in computing the a posteriori error estimator 
% for a decoupled formulation of the quad curl problem
% - Error analysis of a decoupled finite element method for quad-curl problems
%   arXiv:2102.03396
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

Dphi = -3*Dlambda;
ux = u(:,1); uy = u(:,2); uz = u(:,3);
NT = size(elem2face,1);
curlu = zeros(NT,3);
for i = 1:4 % 4 faces
    curlu = curlu ...
       + repmat(ux(elem2face(:,i)),[1,3]).*[zeros(NT,1), Dphi(:,3,i), -Dphi(:,2,i)]...
       + repmat(uy(elem2face(:,i)),[1,3]).*[-Dphi(:,3,i), zeros(NT,1), Dphi(:,1,i)]...
       + repmat(uz(elem2face(:,i)),[1,3]).*[Dphi(:,2,i), -Dphi(:,1,i), zeros(NT,1)];
end