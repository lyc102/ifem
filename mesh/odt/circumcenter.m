function [c,r] = circumcenter(node,elem)
%% CIRCUMCENTER computes circumcenters of a triangulation
%
% [c,r] = circumcenter(node,elem) returns circumcenters of triangles in the
% triangulation represented by (node,elem) along with radius r.
%
% Example
%    theta = [-2*pi/3 -pi/3 0 pi/3 2*pi/3 pi]';
%    node = [cos(theta), sin(theta)];
%    node(end+1,:) = [0,0];
%    elem = delaunayn(node);
%    showmesh(node,elem);
%    c = circumcenter(node,elem);
%    hold on; plot(c(:,1),c(:,2),'r*')
%
% The formula is derived from ODT-based mesh smoothing; see <a href =
% "http://math.uci.edu/~chenlong/thesis.html">Chen Thesis</a> page 52,
% formula (3.23). 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*abs(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
x2 = node(:,1).^2 + node(:,2).^2;
w1 = x2(elem(:,2))+x2(elem(:,3));
w2 = x2(elem(:,1))+x2(elem(:,3));
w3 = x2(elem(:,1))+x2(elem(:,2));
fe1 = [w1 w1].*[ve1(:,2) -ve1(:,1)];
fe2 = [w2 w2].*[ve2(:,2) -ve2(:,1)];
fe3 = [w3 w3].*[ve3(:,2) -ve3(:,1)];
c = 0.25*(fe1 + fe2 + fe3)./[area area];
r = sqrt(sum((node(elem(:,1),:) - c).^2,2));
