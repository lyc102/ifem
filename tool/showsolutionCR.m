function showsolutionCR(node,elem,u,varargin)
%% SHOWSOLUTIONCR plots a CR function u on a triangular mesh in 2-D.
%
%     f = inline('sin(2*pi*x(:,1)).*cos(2*pi*x(:,2))');
%     [node, elem] = squaremesh([0,1,0,1],1/2^4);
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
%% Generate a big triangulation
elemnew = reshape(1:3*NT,NT,3);
nodenew = node(elem(:),:);

%% Evaluate piecewise linear function
elem2edge = dofedge(elem);
unew(1:NT) = -u(elem2edge(:,1)) + u(elem2edge(:,2)) + u(elem2edge(:,3));
unew(NT+(1:NT)) = u(elem2edge(:,1)) - u(elem2edge(:,2)) + u(elem2edge(:,3));
unew(2*NT+(1:NT)) = u(elem2edge(:,1)) + u(elem2edge(:,2)) - u(elem2edge(:,3));

%% Plot the solution
showsolution(nodenew,elemnew,unew,varargin{:});