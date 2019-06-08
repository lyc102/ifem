function [u,A] = surfacePoisson(node,elem,bdEdge,f,g_D,g_N)
%% SURFACEPOISSON Poisson on surfaces.
%
% u = surfacePoisson(node,elem,bdEdge,f,g_D,g_N) assembesl the matrix
% equation Au=b for the linear finite element discritization of Poisson
% equation and on a surface and solves it by multigrid solver. The domain
% is discretized by a triangulation represented by |node| and |elem| and
% the boundary condition is given by |bdEdge|. The function handels |f,
% g_D, g_N| is the data for the surface Poisson equation.
%
% $$-\nabla\cdot(d\nabla u) = f  \quad \Omega$$ 
%
% $$u = g_D \quad \Gamma _D$$
%
% $$\nabla u\cdot n = g_N \quad \Gamma _N$$
%
% [u,A] = surfacePoisson(node,elem,bdEdge,f,g_D,g_N) returns the stiffness
% matrix A which can be used to compute the energy norm as sqrt(u'*A*u).
%
% For a closed surface without boundary, use
% u=Poisson(node,elem,[],@f,[],[])
%
% The only difference between surfacePoisson and Poisson is the way to
% compute the local stiffness matrix.
% 
% Example
%     node = [1,0,0; 0,1,0; -1,0,0; 0,-1,0; 0,0,1; 0,0,-1];
%     elem = [6,1,2; 6,2,3; 6,3,4; 6,4,1; 5,1,4; 5,3,4; 5,3,2; 5,2,1];
%     for i = 1:3
%         [node,elem] = uniformrefine(node,elem);
%     end
%     r = sqrt(node(:,1).^2 + node(:,2).^2 + node(:,3).^2);
%     node = node./[r r r];
%     f = inline('2*p(:,1)','p');
%     [u,A] = surfacePoisson(node,elem,[],f,[],[]);
%     set(gcf,'Units','normal'); set(gcf,'Position',[0,0,0.75,0.45]);
%     subplot(1,2,1);
%     showmesh(node,elem,[130,28],0.65); pause(0.1);
%     subplot(1,2,2); 
%     showsolution(node,elem,u,[130,28]); colorbar;
%
% See also: Poisson, Poisson3
%
% <a href="matlab:ifemdoc surfacePoisson">iFEMdoc surfacePoisson</a>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(node,1); NT = size(elem,1); 
u = zeros(N,1); theta = zeros(NT,3);
%% Compute geometric quantities
ve1 = node(elem(:,3),:)-node(elem(:,2),:);
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
le1 = sqrt(sum(ve1.^2,2));
le2 = sqrt(sum(ve2.^2,2));
le3 = sqrt(sum(ve3.^2,2));
theta(:,1) = pi - acos(sum(ve2.*ve3,2)./(le2.*le3)); 
theta(:,2) = pi - acos(sum(ve3.*ve1,2)./(le3.*le1)); 
theta(:,3) = pi - acos(sum(ve1.*ve2,2)./(le1.*le2)); 
area = 0.5*le1.*le2.*sin(theta(:,3));
clear ve1 ve2 ve3 le1 le2 le3
%% Compute local stiffness matrix
%
% $$A_{ij} = \int _{\tau} \nabla \phi _i\cdot \nabla \phi _j = -0.5\cot(\theta_k),$$
%
% where i,j,k are indices of a flat triangle
%
At = zeros(NT,3,3);
At(:,1,2) = -0.5*cot(theta(:,3));
At(:,1,3) = -0.5*cot(theta(:,2));
At(:,2,3) = -0.5*cot(theta(:,1));
At(:,2,1) = At(:,1,2);
At(:,3,1) = At(:,1,3);
At(:,3,2) = At(:,2,3);
At(:,1,1) = -At(:,1,2)-At(:,1,3);
At(:,2,2) = -At(:,2,1)-At(:,2,3);
At(:,3,3) = -At(:,3,1)-At(:,3,2);
clear theta
%% Assemble stiffness matrix
i = zeros(9*NT,1); j = zeros(9*NT,1); s = zeros(9*NT,1); 
index = 0;
for ti = 1:3
    for tj = 1:3
        i(index+1:index+NT) = elem(:,ti);
        j(index+1:index+NT) = elem(:,tj);
        s(index+1:index+NT) = At(:,ti,tj);
        index = index + NT;
    end
end
A = sparse(i,j,s,N,N);
clear i j s At

%% Assemble right hand side by 3-points quadrature rule
mid1 = (node(elem(:,2),:)+node(elem(:,3),:))/2;
mid2 = (node(elem(:,3),:)+node(elem(:,1),:))/2;
mid3 = (node(elem(:,1),:)+node(elem(:,2),:))/2;
bt1 = area.*(f(mid2)+f(mid3))/6;
bt2 = area.*(f(mid3)+f(mid1))/6;
bt3 = area.*(f(mid1)+f(mid2))/6;
b = accumarray(elem(:),[bt1;bt2;bt3],[N 1]);
clear mid1 mid2 mid3 bt1 bt2 bt3

%% Boundary Conditions
% Find boundary edges and nodes
% case: Poisson(node,elem,[],@f,[],@g_N)
%       Poisson(node,elem,[],@f,@g_D,[])
if isempty(bdEdge) && (~isempty(g_D) || ~isempty(g_N))
    [Dirichlet,Neumann] = findboundary(elem);
end
% case: Poisson(node,elem,bdEdge,@f,...) 
if ~isempty(bdEdge)
	allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
    Dirichlet = allEdge((bdEdge(:) == 1),:);
    Neumann = allEdge((bdEdge(:) == 2),:); 
end
% Dirichlet boundary conditions
if ~isempty(g_D)
    isBdNode = false(N,1); 
    isBdNode(Dirichlet) = true;
    bdNode = find(isBdNode);
    freeNode = find(~isBdNode);
    u(bdNode) = g_D(node(bdNode,:));
    b = b - A*u;
    b(bdNode) = u(bdNode);
end
% Neumann boundary conditions
if (~isempty(g_N) && ~isempty(Neumann))
    freeNode = 1:N;
    ve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
    edgeLength = sqrt(sum(ve.^2,2)); 
    mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
    ge = edgeLength.*(2/3*g_N(mid) + 1/6*g_N(node(Neumann(:,1),:)) ...
                     + 1/6*g_N(node(Neumann(:,1),:)));
    b = b + accumarray([Neumann(:),ones(2*size(Neumann,1),1)], ... 
                   repmat(ge,2,1),[N,1]); 
end
% Pure Neumann boundary condition
if isempty(g_D)
    b = b - mean(b);   % compatilbe condition
    freeNode = 2:N;
    bdNode = 1;
end
clear allEdge Dirichlet Neumann

%% Solve the linear system of algebraic equations

% Multigrid Method for large system
if  N > 4e3
    bdidx = zeros(N,1); 
    bdidx(bdNode) = 1;
    Tbd = sparse(1:N,1:N,bdidx,N,N);
    T = sparse(1:N,1:N,1-bdidx,N,N);
    AD = T*A*T + Tbd;
    u = mg(AD,b,elem);
%%
% Build Dirichlet boundary condition into the matrix A by changing
% |A(bdNode,bdNode)= I|.
end
% Direct solver for small size problems
if ~isempty(freeNode)    
    u(freeNode) = A(freeNode,freeNode)\b(freeNode); 
end
%% TODO: update follow the style of Poisson
