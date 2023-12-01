function h = showsolutionCR(node,elem,u,varargin)
%% SHOWSOLUTIONCR plots a CR function u on a triangular mesh in 2-D.
%
%  Example: a basis 
%
%     [node,elem] = squaremesh([-1,1,-1,1],0.5);
%     edge = myunique(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
%     phi = zeros(size(edge,1),1);
%     phi(6) = 1;
%     showmesh(node,elem,'facecolor','none'); hold on;
%     showsolutionCR(node,elem,phi,[30,26],'facecolor','g','facealpha',0.5,'edgecolor','k');
%
%  Example: quadratic function
%
%     [node,elem] = squaremesh([-1,1,-1,1],0.5);
%     edge = myunique(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
%     f = inline('3 - x(:,1).^2 - x(:,2).^2');
%     u = Lagrangeinterpolate(f,node,elem,'CR',edge);
%     showmesh(node,elem); hold on;
%     showsolutionCR(node,elem,u,[30,26],'facecolor','g','facealpha',0.5,'edgecolor','k');
%
%  Example: sin(pi*x)cos(pi*y)
%
%     f = inline('sin(pi*x(:,1)).*cos(pi*x(:,2))');
%     [node, elem] = squaremesh([0,1,0,1],1/2^2);
%     edge = myunique(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
%     u = Lagrangeinterpolate(f,node,elem,'CR',edge);
%     subplot(1,2,1);
%     showsolutionCR(node,elem,u);
%     subplot(1,2,2);
%     showsolution(node,elem,u,'EdgeColor','k');
%
%  See also showsolution
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

NT = size(elem,1);
%% Generate a big and discontinuous triangulation
elemnew = reshape(1:3*NT,NT,3);
nodenew = node(elem(:),:);

%% Evaluate piecewise linear function
elem2edge = dofedge(elem);
unew(1:NT) = -u(elem2edge(:,1)) + u(elem2edge(:,2)) + u(elem2edge(:,3));
unew(NT+(1:NT)) = u(elem2edge(:,1)) - u(elem2edge(:,2)) + u(elem2edge(:,3));
unew(2*NT+(1:NT)) = u(elem2edge(:,1)) + u(elem2edge(:,2)) - u(elem2edge(:,3));

%% Plot the solution
h = showsolution(nodenew,elemnew,unew,varargin{:});