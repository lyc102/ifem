function [Dlambda,volume,elemSign] = gradbasis3(node,elem)
%% GRADBASIS3 gradient of barycentric basis in 3-D.
%
% [Dlambda,volume,elemSign] = gradbasis3(node,elem) compute gradient of
% barycentric basis and volume of tetrahedron. The array volume is NT by 1
% and volume(t) is the volume of the t-th tetrahedron. Dlambda(1:NT, 1:3,
% 1:4) the first index is the label of tetrahedron, the second is the x-y-z
% coordinate, and the last one is the label of four indices of a tet. For
% example, Dlambda(t,:,1) is the gradient of lambda of the 1st index of the
% t-th tet. The elemSign array taking values 1 or -1 records the
% sign of the signed volume.
% 
% See also gradbasis, gradu, gradu3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

NT = size(elem,1);
face = uint32([elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])]);
v12 = node(face(:,2),:)-node(face(:,1),:);
v13 = node(face(:,3),:)-node(face(:,1),:);
normal = mycross(v12,v13,2);
v12 = v12(3*NT+1:4*NT,:); 
v13 = v13(3*NT+1:4*NT,:);
v14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = dot(mycross(v12,v13,2),v14,2)/6;
Dlambda = zeros(NT,3,4);
Dlambda(1:NT,:,1) = normal(1:NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda(1:NT,:,2) = normal(NT+1:2*NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda(1:NT,:,3) = normal(2*NT+1:3*NT,:)./[6*volume, 6*volume, 6*volume];
Dlambda(1:NT,:,4) = normal(3*NT+1:4*NT,:)./[6*volume, 6*volume, 6*volume];

% When the tetrahedron is not positive orientated, we reverse the sign of
% the volume. The sign of Dlambda is always right since signed volume is
% used in the computation. 
idx = (volume<0); 
volume(idx,:) = -volume(idx,:);
elemSign = ones(NT,1);
elemSign(idx) = -1;