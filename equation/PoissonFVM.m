function [u,AD,F,N,itStep,cputime] = PoissonFVM(node,elem,pde,bdEdge,centerpoint,solver)
%% POISSONFVM Poisson equation: linear finite volume method.
%
%   u = PoissonFVM(node,elem,pde,bdEdge,centerpoint) produces the linear
%   finite element approximation of the Poisson equation
%
%       -div(d*grad(u))=f  in \Omega, with
%       Dirichlet boundary condition u=g_D on \Gamma_D,
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
%   The mesh is given by node and elem and the boundary edge is given by
%   bdEdge. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, g_R, or d.
%   For general elliptic equations with convection and reaction
%   coefficients, see ellipticpde.
%
%   The function PoissonFVM assembes the matrix equation AD*u = F and solves
%   it by the direct solver (small size < 2e3) or the multigrid solver
%   (large size >=2e3). The Dirichlet boundary condition is built into the
%   matrix AD and the Neumann boundary condition is build into F.
%
%   The diffusion coefficient d is a scalar function or a column array with
%   the same length as elem array.
%
%   The type of linear finite volume is given by centerpoint, which can be
%   'barycenter' and 'longedgecenter', corresponds to the following two
%   types of choices of control volumes, respectively.
%   • Type A: c_t is the barycenter of t                 --- barycenter
%   • Type B: c_t is the middle point of the longest edge--- longedgecenter
%
%   [u,AD,F] = PoissonFVM(node,elem,pde,bdEdge, centerpoint) returns the
%   modified stiffness matrix AD and the load vector F. The
%   solution u = AD\F. The matrix AD can be used to test other solvers.
%
%   [u,AD,F,A] = PoissonFVM(node,elem,pde,bdEdge) returns also the
%   non-modified stiffness matrix A,  The kernel of A consists of constant
%   vectors.
%
%   u = PoissonFVM(node,elem,pde,bdEdge,solver) specifies the solver by the
%   string 'solver'. If solver == 'direct', the Matlab built in direct
%   solver \ (mldivide) is used. If solver == 'mg', multigrid-type solvers
%   mg is used. If solver == 'notsolve', the solution u = []. The default
%   setting is to use the direct solver for small size problems (< 2e3) and
%   multigrid solvers for large size problems.
%
%   Example
%
%         squarePoissonFVM
%
%   See also Poisson Poisson3, squarePoisson, Lshape, crack, mg
%
%   Created by Ming Wang, and merged two types FVM at May 15, 2011.
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(node,1); NT = size(elem,1);

%% Diffusion coefficient
if isfield(pde,'d') && ~isempty(pde.d)
   if isnumeric(pde.d)
      K = pde.d;                      % d is an array
   else                            % d is a function
      center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
      K = pde.d(center);  
   end
else
    K=[];
end

%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
A = sparse(N,N);
for i = 1:3
    for j = i:3
        %%
        % $A_{ij}|_{\tau} = \int_{\tau}K\nabla \phi_i\cdot \nabla \phi_j
        % dxdy$ 
        Aij = (Dphi(:,1,i).*Dphi(:,1,j)+Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if isfield(pde,'d') && ~isempty(pde.d)
          Aij = K.*Aij;
        end
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Aij; Aij],N,N);        
        end        
    end
end
clear  Aij

%% Assemble the right hand side
F = zeros(N,1);
if isfield(pde,'f') && isreal(pde.f) && (pde.f==0)
    pde.f = [];
end
if isfield(pde,'f') && ~isempty(pde.f)
    mid =zeros(NT,2,3);
    mid(:,:,1) = (node(elem(:,2),:)+node(elem(:,3),:))/2;
    mid(:,:,2) = (node(elem(:,3),:)+node(elem(:,1),:))/2;
    mid(:,:,3) = (node(elem(:,1),:)+node(elem(:,2),:))/2;
    ft=zeros(NT,3); locId = [ 2 3; 3 1; 1 2];
    [lambda,weight] = quadpts(2); nQuad = size(lambda,1);
    switch (centerpoint)
        case 'barycenter'
            center = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
            for p = 1:nQuad
                % quadrature points in the x-y coordinate
                for i=1:3
                    pxy1 = lambda(p,1)*node(elem(:,i),:) ...
                        + lambda(p,2)*mid(:,:,locId(i,1)) ...
                        + lambda(p,3)*center;
                    pxy2 = lambda(p,1)*node(elem(:,i),:) ...
                        + lambda(p,2)*mid(:,:,locId(i,2)) ...
                        + lambda(p,3)*center;
                    ft(:,i)=ft(:,i)+weight(p)*(pde.f(pxy1)+pde.f(pxy2)).*area/6;
                end
            end
        case 'longedgecenter'
            for p = 1:nQuad
                % quadrature points in the x-y coordinate
                for i = 1:3
                    pxy = lambda(p,1)*node(elem(:,i),:) ...
                        + lambda(p,2)*mid(:,:,locId(i,1)) ...
                        + lambda(p,3)*mid(:,:,locId(i,2));
                    ft(:,i) = ft(:,i)+weight(p)*pde.f(pxy).*area/4;
                end
                pxy = lambda(p,1)*mid(:,:,1) ...
                    + lambda(p,2)*mid(:,:,2) ...
                    + lambda(p,3)*mid(:,:,3);
                ft(:,1) = ft(:,1)+weight(p)*pde.f(pxy).*area/4;
            end
    end
    F = accumarray(elem(:),ft(:),[N 1]);
end

%% Boundary Conditions
if nargin<=3, bdEdge = []; end
[AD,F,u,freeNode] = getbdFVM(A,F);

%% Solve the system of linear equations
if (nargin<=5)
    if (N < 2e3)    % Direct solver for small size systems
        solver = 'direct';
    else            % Multigrid-type solver for large size systems
        solver = 'mg';
    end
end

if strcmp(solver,'direct') && ~isempty(freeNode)
    u(freeNode) = AD(freeNode,freeNode)\F(freeNode);         
elseif strcmp(solver,'mg')
%     option.tol=1e-12; option.printlevel = 0;
option.tol=1e-12;
     [u,tempvar,itStep,tempvar,cputime]  = mg(AD,F,elem,option);%#ok<*ASGLU>
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,F,u,freeNode] = getbdFVM(A,F)
    %% GETBDFVM Boundary conditions for Poisson equation: P1 linear FVM.
    %
    % Reference: Long Chen. Lecture Notes.
    %
    % Copyright (C) Long Chen. See COPYRIGHT.txt for details.

    N = size(node,1); 
    u = zeros(N,1); 

    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end
    if (isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R))
        bdEdge = [];
    end

    %% Find boundary nodes and edges
    Dirichlet = []; Neumann = []; Robin = [];
    if isempty(bdEdge) % case: bdEdge = [], one bd condition is non-empty
        if (~isempty(pde.g_D) || ~isempty(pde.g_N) || ~isempty(pde.g_R))
            % one type of boundary condition
            [Dirichlet,Neumann] = findboundary(elem);
            if ~isempty(pde.g_R)
                Robin = Neumann;
            end
        end
        % case: bdEdge = [], pde.g_D = pde.g_N = pde.g_R = [];
        % it is equivalent to homogenous Neumann boundary condition
    else  % case: bdEdge gives bd condition
        % mixed boundary conditions
        allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        Dirichlet = allEdge((bdEdge(:) == 1),:);
        Neumann = allEdge((bdEdge(:) == 2) | (bdEdge(:) == 3),:);
        Robin = allEdge((bdEdge(:) == 3),:);
    end

    isBdNode = false(N,1); 
    isBdNode(Dirichlet(:)) = true;
    bdNode = find(isBdNode);
    freeNode = find(~isBdNode);

    %% Neumann boundary condition
    if ~isempty(pde.g_N) && isnumeric(pde.g_N) && (pde.g_N==0)
        pde.g_N = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.g_N)
        ve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
        edgeLen = sqrt(sum(ve.^2,2)); 
        mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
        nQuad = 4;[lambda,weight] = quadpts1(nQuad);
        ge = zeros(size(Neumann,1),2);
        for ip = 1:nQuad
          pxy=lambda(ip,1)*node(Neumann(:,1),:)+lambda(ip,2)*mid;               
          ge(:,1)=ge(:,1)+weight(ip)*edgeLen.*pde.g_N(pxy)/2;
          pxy=lambda(ip,1)*node(Neumann(:,2),:)+lambda(ip,2)*mid;               
          ge(:,2)=ge(:,2)+weight(ip)*edgeLen.*pde.g_N(pxy)/2;
        end
        F = F + accumarray(Neumann(:), ge(:),[N,1]); 
    end
    % The case: non-empty Neumann but g_N=0 or g_N=[] corresponds to the zero
    % flux boundary condition on Neumann edges and no modification is needed.

    %% Dirichlet boundary condition
    if ~isempty(pde.g_D) && isnumeric(pde.g_D) && (pde.g_D==0)
        pde.g_D = [];
    end
    if ~isempty(Dirichlet) && ~isempty(pde.g_D)
        u(bdNode) = pde.g_D(node(bdNode,:));
        F = F - A*u;
        F(bdNode) = u(bdNode);
    end
    % The case: non-empty Dirichlet but g_D=0 or g_D=[] corresponds to the zero
    % Dirichlet boundary condition and no modification of u,F is needed.

    %% Robin boundary condition
    if ~isempty(pde.g_R) && isnumeric(pde.g_R) && (pde.g_R==0)
        pde.g_R = [];
    end
    if ~isempty(Robin) && ~isempty(pde.g_R)
        ve = node(Robin(:,1),:) - node(Robin(:,2),:);
        edgeLength = sqrt(sum(ve.^2,2)); 
        mid = (node(Robin(:,1),:) + node(Robin(:,2),:))/2;
        ii = [Robin(:,1),Robin(:,1),Robin(:,2),Robin(:,2)];
        jj = [Robin(:,1),Robin(:,2),Robin(:,1),Robin(:,2)];
        temp = pde.g_R(mid).*edgeLength;
        ss = [1/3*temp, 1/6*temp, 1/6*temp, 1/3*temp];
        A = A + sparse(ii,jj,ss,N,N);
    end

    %% Pure Neumann boundary condition
    if isempty(Dirichlet) && isempty(Robin)
        F = F - mean(F);   % compatilbe condition: sum(F) = 0
        freeNode = 2:N;    % eliminate the kernel by enforcing u(1) = 0;
        bdNode = 1;
    end

    %% Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % |AD(bdNode,bdNode)=I, AD(bdNode,freeNode)=0, AD(freeNode,bdNode)=0|.
    bdidx = zeros(N,1); 
    bdidx(bdNode) = 1;
    Tbd = spdiags(bdidx,0,N,N);
    T = spdiags(1-bdidx,0,N,N);
    AD = T*A*T + Tbd;
    %% TODO - check the logic of bd condition - write help
    end
end
