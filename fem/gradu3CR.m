function [Du,volume,Dphi] = gradu3CR(node,elem,elem2face,u,Dlambda)
%% GRADU3CR gradient of a Crouzeix-Raviart finite element function.
%
% Du = GRADU3CR(node,elem,elem2face,u) compute the gradient of C-R
% u on a mesh representing by (node,elem).
% 
% [Du,volume,Dphi] = GRADUCR(node,elem,elem2face,u) also outputs area and
% Dphi which is the gradient of nonconforming P1 basis (1-3\lambda). 
%
% Du = GRADUCR(node,elem,elem2face,u,Dlambda) compute the gradient with Dlambda. 
%
% See also gradu3, gradbasis3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('Dlambda','var')
    [Dlambda,volume] = gradbasis3(node,elem);
end
Dphi = -3*Dlambda;
Du = dot(permute(repmat(u(elem2face),[1,1,3]),[1,3,2]),Dphi,3);