function [Du,area,Dlambda] = graduCR(node,elem,elem2edge,u,Dlambda)
%% GRADUCR gradient of a finite element function.
%
% Du = GRADUCR(node,elem,elem2edge,u) compute the gradient of a finite element function
% u on a mesh representing by (node,elem).
% 
% [Du,area,Dlambda] = GRADUCR(node,elem,elem2edge,u) also outputs area and Dlambda
% which is the gradient of P1 basis. 
%
% Du = GRADUCR(node,elem,u,Dlambda) compute the gradient with Dlambda. It
% will save time when Dlambda is available. See recovery.m
%
% See also gradu3, gradbasis
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('Dlambda','var')
    [Dlambda,area] = gradbasis(node,elem);
end
dudx =  u(elem2edge(:,1)).*Dlambda(:,1,1) + u(elem2edge(:,2)).*Dlambda(:,1,2) ...
      + u(elem2edge(:,3)).*Dlambda(:,1,3);
dudy =  u(elem2edge(:,1)).*Dlambda(:,2,1) + u(elem2edge(:,2)).*Dlambda(:,2,2) ...
      + u(elem2edge(:,3)).*Dlambda(:,2,3);         
Du = -2*[dudx, dudy];