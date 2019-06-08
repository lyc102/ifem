function [Clambda,area,elemSign] = curlbasis(node,elem)
%% CURLBASIS gradient of barycentric basis. 
%
% [Clambda,area,elemSign] = CURLBASIS(node,elem) compute gradient of
% barycentric basis and areas of triangles. The array area is NT by 1 and
% area(t) is the volume of the t-th tetrahedron. Clambda(1:NT, 1:3, 1:3)
% the first index is the label of tetrahedron, the second is the x-y
% coordinate, and the last one is the label of three indices of a triangle.
% For example, Clambda(t,:,1) is the gradient of lambda of the 1st index of
% the t-th triangle. The elemSign array taking values 1 or -1 records the
% sign of the signed area.
%
% See also gradbasis, gradu, gradu3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

NT = size(elem,1);
% $\nabla \phi_i = rotation(l_i)/(2|\tau|)$
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
Clambda(1:NT,:,1) = [ve1(:,1)./(2*area), ve1(:,2)./(2*area)];
Clambda(1:NT,:,2) = [ve2(:,1)./(2*area), ve2(:,2)./(2*area)];
Clambda(1:NT,:,3) = [ve3(:,1)./(2*area), ve3(:,2)./(2*area)];

% When the triangle is not positive orientated, we reverse the sign of the
% area. The sign of Dlambda is always right since signed area is used in
% the computation.
idx = (area<0); 
area(idx,:) = -area(idx,:);
elemSign = ones(NT,1);
elemSign(idx) = -1;