function [u,elem2edge,A] = PoissonCR(node,elem,bdEdge,f,g_D,g_N)
% POISSONCR assembel the matrix Ax=b for the Crouzeix and Raviart nonconforming
% linear finite element discritization of Poisson equation and solve it by 
% direct solver
%          $-\Delta u = f $    in  $\Omega$ 
%             $u = g_D$ on  $\Gamma _D$ Dirichelet boundary  
%         $\nabla u\cdot n = g_N$ on $\Gamma _N$ Neumann boundary 
% in a domain described by node and elem, with boundary conditions 
%
%    u = Poisson(node,elem,bdEdge,f,g_D,g_N)  % mixed boundary condition
%   [u,A,b] = Poisson(node,elem,bdEdge,f,g_D,g_N) 

% Input
%    bdEdge represents the boundary edges. It can be omitted for 
%    u = PoissonCR(node,elem,[],f,g_D,[])       % Dirichlet only 
%    u = PoissonCR(node,elem,[],f,[],g_N)       % Neumann only
%    u = PoissonCR(node,elem,[],f,[],[])        % zero flux
%
%    f, g_D, g_N: the load and boundary conditions

%--------------------------------------------------------------------------
% Copyright Long Chen.
%--------------------------------------------------------------------------

%-------- Construct data structure ---------------------------------------- 
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[edge, i2, jj] = unique(totalEdge,'rows');
elem2edge = reshape(jj,size(elem,1),3);
%-------- Construct A and u ----------------------------------------
NT = size(elem,1); NE = size(edge,1); A = sparse(NE,NE); u = zeros(NE,1);
%-------- Compute vedge, edge as a vector, and area of each element ------
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
%-------- Assemble stiffness matrix ---------------------------------------
for i = 1:3
    for j = i:3
        Aij = (ve(:,1,i).*ve(:,1,j) + ve(:,2,i).*ve(:,2,j))./area;
        if (j==i)
            A = A + sparse(elem2edge(:,i),elem2edge(:,j),Aij,NE,NE);
        else
            A = A + sparse([elem2edge(:,i);elem2edge(:,j)], ...
                         [elem2edge(:,j);elem2edge(:,i)],[Aij; Aij],NE,NE);        
        end        
    end
end
clear ve
%-------- Assembing right hand side by 3-points rule ----------------------
mid1 = (node(elem(:,2),:) + node(elem(:,3),:))/2;
mid2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;
mid3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;
bt1 = area.*f(mid1)/3;
bt2 = area.*f(mid2)/3;
bt3 = area.*f(mid3)/3;
b = accumarray(elem2edge(:),[bt1;bt2;bt3],[NE 1]);
clear mid1 mid2 mid3 bt1 bt2 bt3
%-------- Set boundary edges and nodes ------------------------------------
if ~isempty(bdEdge)
    Dirichlet = elem2edge(bdEdge(:)==1);
    isNeumann = elem2edge(bdEdge(:)==2);
    Neumann = edge(isNeumann,:);
else
    i1(jj(3*NT:-1:1)) = 3*NT:-1:1; 
    i1=i1'; 
    Dirichlet = find(i1 == i2);
    Neumann = edge((i1 == i2),:);
    isNeumann = (i1==i2);
end
%-------------------- Neumann boundary conditions -------------------------
if (~isempty(g_N) && ~isempty(Neumann))
    Nve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
    edgeLength = sqrt(sum(Nve.^2,2)); 
    mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
    b(isNeumann) = b(isNeumann) + edgeLength.*g_N(mid);
end
%-------------------- Dirichlet boundary conditions------------------------
if ~isempty(g_D)
    isBdEdge = false(NE,1); 
    isBdEdge(Dirichlet) = true;
    freeEdge = find(~isBdEdge);
    mid = (node(edge(Dirichlet,1),:) + node(edge(Dirichlet,2),:))/2; 
    u(Dirichlet) = g_D(mid);
    b = b - A*u;
else
    b = b - mean(b); % compatible condition for pure Neumann bd condition
    freeEdge = 2:NE;  % pure Neumann boundary condition has a kernel
end
%-------- Solve the linear system Au = b for the free edges ------------ 
if size(freeEdge)>0
    u(freeEdge) = A(freeEdge,freeEdge)\b(freeEdge);
end