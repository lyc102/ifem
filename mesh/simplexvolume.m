function [v,elemSign] = simplexvolume(node,elem)
%% SIMPLEXVOLUME Simplex volume
%
%   v = SIMPLEXVOLUME(node,elem) computes signed measure of simplicies in the
%   triangulation given by node and elem. 
%
%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
%   Modified by Long Chen Feb 05, 2010.

NT = size(elem,1);
switch size(elem,2)-1 % dimension of simplex
 case 1  % an interval
    d12 = node(elem(:,2),:)-node(elem(:,1),:);
    v = d12;
 case 2  % a triangle
    d12 = node(elem(:,2),:)-node(elem(:,1),:);
    d13 = node(elem(:,3),:)-node(elem(:,1),:);
    v = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;
    if size(node,2) == 3     % surface triangles
        normal = mycross(d12,d13,2);
        v = sqrt(sum(normal.^2,2))/2;
    end
 case 3  % a tetrahedron
    d12 = node(elem(:,2),:)-node(elem(:,1),:);
    d13 = node(elem(:,3),:)-node(elem(:,1),:);
    d14 = node(elem(:,4),:)-node(elem(:,1),:);
    v = dot(mycross(d12,d13,2),d14,2)/6;
end

% When the tetrahedron is not positive orientated, we reverse the sign of
% the volume. The sign of Dlambda is always right since signed volume is
% used in the computation. 
idx = (v<0); 
v(idx,:) = -v(idx,:);
elemSign = ones(NT,1);
elemSign(idx) = -1;