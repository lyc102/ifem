function [Du,volume,Dlambda] = gradu3(node,elem,u,Dlambda)
%% GRADU3 gradient of a finite element function in 3-D.
%
% Du = gradu3(node,elem,u) compute the gradient of a finite element function
% u on a 3-D mesh representing by (node,elem).
% 
% [Du,volume,Dlambda] = gradu3(node,elem,u) also outputs area and Dlambda
% which is the gradient of P1 basis. 
%
% Du = gradu3(node,elem,u,Dlambda) compute the gradient with Dlambda. It
% will save time when Dlambda is available. See recovery3.m
%
% See also gradu3, gradbasis
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if (nargin <= 3)
    [Dlambda,volume] = gradbasis3(node,elem);
end
dudx = u(elem(:,1)).*Dlambda(:,1,1)+u(elem(:,2)).*Dlambda(:,1,2) ...
      +u(elem(:,3)).*Dlambda(:,1,3)+u(elem(:,4)).*Dlambda(:,1,4);
dudy = u(elem(:,1)).*Dlambda(:,2,1)+u(elem(:,2)).*Dlambda(:,2,2) ...
      +u(elem(:,3)).*Dlambda(:,2,3)+u(elem(:,4)).*Dlambda(:,2,4);
dudz = u(elem(:,1)).*Dlambda(:,3,1)+u(elem(:,2)).*Dlambda(:,3,2) ...
      +u(elem(:,3)).*Dlambda(:,3,3)+u(elem(:,4)).*Dlambda(:,3,4);
Du = [dudx, dudy, dudz];