function [u,A,b] = PoissonS(node,elem,Dirichlet,Neumann,f,g_D,g_N)
%% POISSONS solve the 2-D Poisson equation by linear finite element
%
%     -\Delta u = f,
%             u = g_D on the Dirichelet boundary edges
%         du/dn = g_N on the Neumann boundary edges
% in a domain described by node and elem, with boundary edges Dirichlet,
% Neumann. It is a simplified version of Poisson.
%
% Input: 
%   node, elem: standard mesh data;
%   Dirichlet, Neumann: boundary edges;
%   f: functional handle right side or data
%   g_D: functional handle for Dirichelet condition
%   g_N: functional handle for Neumann condition
%
% Output:
%   u: solution on the current mesh
%   A: stiffness matrix
%
% Example
%
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     f = inline('2*pi^2*cos(pi*p(:,1)).*cos(pi*p(:,2))','p');
%     exactu = inline('cos(pi*p(:,1)).*cos(pi*p(:,2)) - 1','p');
%     err = zeros(4,1); N = zeros(4,1);
%     for k = 1:4
%       [node,elem] = uniformrefine(node,elem);
%       Dirichlet = findboundary(elem);
%       [u,A] = PoissonS(node,elem,Dirichlet,[],f,exactu);
%       uI = exactu(node);
%       err(k) = sqrt((u-uI)'*A*(u-uI));
%       N(k) = length(u);
%       figure(1); showresult(node,elem,u);
%     end
%     figure; showrate(N,err);     
%
% See also Poisson, Poisson3
%
% The code is based on the following reference but optimized using
% vectorization to avoid for loops.
%
% Reference:
%
%    Jochen Alberty, Carsten Carstensen, Stefan Funken, Remarks Around 50
%    Lines of MATLAB: Short Finite Element Implementation, Numerical
%    Algorithms, Volume 20, pages 117-137, 1999.
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(node,1);
A = sparse(N,N); u = zeros(N,1);

%% Assembing stiffness matrix
ve(:,:,3) = node(elem(:,2),:) - node(elem(:,1),:);
ve(:,:,1) = node(elem(:,3),:) - node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:) - node(elem(:,3),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
for i = 1:3
    for j = 1:3
        Aij = (ve(:,1,i).*ve(:,1,j) + ve(:,2,i).*ve(:,2,j))./(4*area);
        A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
    end
end

%% Assembing right side
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
b = accumarray(elem(:),repmat(f(center).*area/3,3,1),[N,1]);

%% Neumann boundary conditions
if(~isempty(Neumann))
    d = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2)); 
    mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
    b = b + accumarray(Neumann(:),repmat(d.*g_N(mid)/2,2,1),[N,1]);
end

%% Dirichlet boundary conditions
isBdNode = false(N,1); 
if ~isempty(Dirichlet)
    isBdNode(Dirichlet(:)) = true;
    u(isBdNode) = g_D(node(isBdNode,:));
    b = b - A*u;
else % Pure Neumann boundary condition
    b = b - mean(b);
    isBdNode(1) = true;
end

%% Solve Au = b
freeNode = find(~isBdNode);
u(freeNode) = A(freeNode,freeNode)\b(freeNode);