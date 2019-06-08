function [u,edge,A,F,Ndof,itStep,cputime] = PoissonP2FVM(node,elem,bdEdge,pde,centerpoint,solver,varargin)
%% POISSONP2FVM Poisson equation: P2 quadratic finite volume method.
%
%  u = PoissonP2FVM(node,elem,bdEdge,f,g_D,g_N,g_R,d) produces a quadratic
%  finite volume approximation of the Poisson equation in two dimensions.
%
%  See also PoissonP2FVMc, PoissonP2, Poisson3P2 
%
%  Modified by Long Chen and Ming Wang. 
%
%  Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Construct Data Structure
[elem2dof,edge,bdDof] = dofP2(elem);
N = size(node,1);  NT = size(elem,1);  Ndof = N+size(edge,1);

%% compute diffusion coefficient
if isfield(pde,'d')
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
[Dlambda,area] = gradbasis(node,elem);

%% Compute element-wise stiffness matrix
AT = zeros(NT,3,3);
for i = 1:3
    for j = 1:3
            AT(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2).*area;
            if ~isempty(K)
                AT(:,i,j) = K.*AT(:,i,j);
            end
    end
end
clear Dlambda
BT = -4/3*AT;
DT = 8/3*AT;
for i=1:3
	DT(:,i,i) = 4/3*(AT(:,1,1) + AT(:,2,2) + AT(:,3,3));
end
c1 = -2*AT(:,2,3); c2 = -2*AT(:,1,3); c3 = -2*AT(:,1,2);
switch (centerpoint)
    case 'barrycenter'
        CT(1:NT,1,1:3) = -[-c1/3 + c2/2 + c3/2, c2/6 - c3/2, c3/6 - c2/2];
        CT(1:NT,2,1:3) = -[c1/6 - c3/2, c1/2 - c2/3 + c3/2, c3/6 - c1/2];
        CT(1:NT,3,1:3) = -[c1/6 - c2/2, c2/6 - c1/2, c1/2 + c2/2 - c3/3];
    case 'longedgecenter'
        CT(1:NT,1,1:3) = 1/2*[-c2-c3+2*c1, c2+c3, c2+c3];
        CT(1:NT,2,1:3) = 1/2*[-c1+c3, -c1-c3, c1-c3];
        CT(1:NT,3,1:3) = 1/2*[-c1+c2, c1-c2, -c1-c2];
end

%% Assemble global stiffness matrix
A = sparse(Ndof,Ndof);
for i = 1:3
    for j = 1:3
        elem2dof = double(elem2dof);
        A = A + sparse([elem2dof(:,i); elem2dof(:,i); elem2dof(:,i+3);elem2dof(:,i+3)],...
            [elem2dof(:,j); elem2dof(:,j+3); elem2dof(:,j);elem2dof(:,j+3)],...
            [AT(:,i,j); CT(:,i,j); BT(:,i,j); DT(:,i,j)],Ndof,Ndof);
    end
end
clear AT BT CT DT

%% Assembing right hand side
F = zeros(Ndof,1);
if isfield(pde,'f') && isreal(pde.f) && (pde.f==0)
    pde.f = [];
end
if isfield(pde,'f') && ~isempty(pde.f)
    mid =zeros(NT,2,3);
    mid(:,:,1) = (node(elem(:,2),:)+node(elem(:,3),:))/2;
    mid(:,:,2) = (node(elem(:,3),:)+node(elem(:,1),:))/2;
    mid(:,:,3) = (node(elem(:,1),:)+node(elem(:,2),:))/2;
    ft=zeros(NT,3);fe=zeros(NT,3);locId = [ 2 3; 3 1; 1 2];
    [lambda,weight] = quadpts(2); nQuad = size(lambda,1);
    phi(:,1) = 4*lambda(:,2).*lambda(:,3); % bubble function on middle point
    phi(:,2) = 4*lambda(:,3).*lambda(:,1);
    phi(:,3) = 4*lambda(:,1).*lambda(:,2);
    % quadrature for linear bases part -- different for two kinds of FVM.
    switch (centerpoint)
        case 'barrycenter'
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
                % quadrature for linear element control volume
                for i = 1:3 % associate the integral on subtriangle to each vertex.
                    pxy = lambda(p,1)*node(elem(:,i),:) ...
                        + lambda(p,2)*mid(:,:,locId(i,1)) ...
                        + lambda(p,3)*mid(:,:,locId(i,2));
                    ft(:,i) =ft(:,i)+weight(p)*pde.f(pxy).*area/4;
                end
                pxy = lambda(p,1)*mid(:,:,1) ... % associate the integral on the
                    + lambda(p,2)*mid(:,:,2) ... % middle triangle to the first vertex
                    + lambda(p,3)*mid(:,:,3);
                ft(:,1) = ft(:,1)+weight(p)*pde.f(pxy).*area/4;
            end
    end
    % quadrature for quadratic bubble bases part -- the same 
    for p = 1:nQuad
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        for i = 1:3
            fe(:,i) = fe(:,i) + weight(p)*phi(p,i)*pde.f(pxy).*area;
        end
    end
    % accumarray
    F = accumarray(elem2dof(:),[ft(:);fe(:)],[Ndof 1]);
    
end
%% Boundary Conditions
if nargin<=3, bdEdge = []; end
[AD,F,u,freeNode] = getbdP2FVM(A,F);

%% Solve the system of linear equations
% set solver type
if (nargin<=5)
    if (N < 1e3)    % Direct solver for small size systems
        solver = 'direct';
    else            % Multigrid-type solver for large size systems
        solver = 'mg';
    end
end
% solve
if strcmp(solver,'direct') && ~isempty(freeNode)
    u(freeNode) = AD(freeNode,freeNode)\F(freeNode);
elseif strcmp(solver,'mg')
%     option.solver = 'gmres';
    option.solver = 'Vcycle';
%     [u,tempvar,itStep,tempvar,cputime] = mg(AD,F,elem,option,edge,'HB');%#ok<*ASGLU>
    [u,info] = mg(AD,F,elem,option,'HB');
    itStep = info.itStep; cputime = info.time;
end
% Change bases from HB to NB.
u(N+1:end)=u(N+1:end)+(u(edge(:,1))+u(edge(:,2)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,F,u,freeDof] = getbdP2FVM(A,F)
    %
    %% Set up boundary and basic parameter
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'d'), pde.d = []; end
    N = size(node,1); Ndof = N +size(edge,1);u = zeros(Ndof,1);

    %% Find boundary nodes and edges
    if ~isempty(bdEdge)     
        elem2edge=elem2dof(:,4:6)-N;
        isDirichlet = unique(elem2edge(bdEdge(:)==1));
        isNeumann   = unique(elem2edge((bdEdge(:)==2)|(bdEdge(:) == 3)));
        isRobin     = unique(elem2edge(bdEdge(:)==3));
        Dirichlet = edge(isDirichlet,:);
        Neumann   = edge(isNeumann,:);
        Robin   = edge(isRobin,:);
        Dirichlet = [unique(Dirichlet(:)); N+isDirichlet];
    else % bdEdge = []
        if ~isempty(pde.g_D)     % case: pure Dirichlet
            Dirichlet = bdDof;
            Neumann = [];
            Robin = [];
        elseif ~isempty(pde.g_N) % case: pure Neumann
            Dirichlet = [];
            Neumann = edge(bdDof>N,:);
            Robin = [];
        elseif ~isempty(pde.g_R) % case: Robin BC
            Robin = edge(bdDof>N,:);
        else                     % case: Homogenous Neumann BC
            Dirichlet = [];
            Neumann = [];
            Robin = [];
        end
    end

    isbdDof = false(Ndof,1); 
    isbdDof(Dirichlet) = true;
    bdDof   = Dirichlet;
    freeDof = find(~isbdDof);

    %% Neumann boundary condition
    if (~isempty(pde.g_N) && ~isempty(Neumann))
        nQuad = 4;[lambda,weight] = quadpts1(nQuad);
        % bases and length of edge
        phi= 4*lambda(:,1).*lambda(:,2); % bubble func
        ve = node(Neumann(:,1),:) - node(Neumann(:,2),:);
        edgeLen = sqrt(sum(ve.^2,2));
        % linear part
        mid = (node(Neumann(:,1),:) + node(Neumann(:,2),:))/2;
        ge = zeros(length(isNeumann),2);
        for ip = 1:nQuad
            pxy=lambda(ip,1)*node(Neumann(:,1),:)+lambda(ip,2)*mid;
            ge(:,1)=ge(:,1)+weight(ip)*edgeLen.*pde.g_N(pxy)/2;
            pxy=lambda(ip,1)*node(Neumann(:,2),:)+lambda(ip,2)*mid;
            ge(:,2)=ge(:,2)+weight(ip)*edgeLen.*pde.g_N(pxy)/2;
        end   
        F(1:N) = F(1:N) + accumarray(Neumann(:), ge(:),[N,1]);
        % quadratic part
        for ip = 1:nQuad
          pxy = lambda(ip,1)*node(Neumann(:,1),:)+lambda(ip,2)*node(Neumann(:,2),:);               
          F(N+isNeumann)=F(N+isNeumann)+weight(ip)*phi(ip)*edgeLen.*pde.g_N(pxy);
        end       
    end


    %% Dirichlet boundary conditions
    if ~isempty(pde.g_D) && isnumeric(pde.g_D) && (pde.g_D==0)
        pde.g_D = [];
    end
    if (~isempty(pde.g_D)&& ~isempty(Dirichlet))
        idx = (bdDof > N);                              % index of edge nodes
        u(bdDof(~idx)) = pde.g_D(node(bdDof(~idx),:));  % bd value at vertex dofs
        bdEdgeIdx = bdDof(idx)-N;
        xybdEdgeNode = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2; 
        u(bdDof(idx)) = pde.g_D(xybdEdgeNode)-(u(edge(bdEdgeIdx,1))+u(edge(bdEdgeIdx,2)))/2;
        F = F - A*u;
        F(bdDof) = u(bdDof);
    end


    %% Robin boundary condition
    if ~isempty(Robin) && ~isempty(g_R)
        % TODO: change to FVM bases
        ve = node(Robin(:,1),:) - node(Robin(:,2),:);
        edgeLen = sqrt(sum(ve.^2,2)); 
        mid = (node(Robin(:,1),:) + node(Robin(:,2),:))/2;
        ii = [Robin(:,1),Robin(:,1),Robin(:,2),Robin(:,2)];
        jj = [Robin(:,1),Robin(:,2),Robin(:,1),Robin(:,2)];
        temp = g_R(mid).*edgeLen;
        ss = [1/3*temp, 1/6*temp, 1/6*temp, 1/3*temp];
        A = A + sparse(ii,jj,ss,N,N);
    end

    %% Pure Neumann boundary condition
    if isempty(Dirichlet)
        F = F - mean(F);   % compatilbe condition: sum(F) = 0
        freeDof = 2:Ndof;  % eliminate the kernel by enforcing u(1) = 0;
    end
    %% 
    % *Modify the matrix*
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % |AD(bdDof,bdDof)=I, AD(bdDof,freeDof)=0, AD(freeDof,bdDof)=0|.

    bdidx = zeros(Ndof,1); 
    bdidx(bdDof) = 1;
    Tbd = sparse(1:Ndof,1:Ndof,bdidx,Ndof,Ndof);
    T = sparse(1:Ndof,1:Ndof,1-bdidx,Ndof,Ndof);
    AD = T*A*T + Tbd;
    end
end
