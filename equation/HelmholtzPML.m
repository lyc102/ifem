function [u,eqn,info] = HelmholtzPML(node,elem,pde,bdFlag,option)
%% HELMHOLTZ Poisson equation: P1 linear element.
%
%   u = Helmholtz(node,elem,pde,bdFlag) produces the linear finite element
%   approximation of the Poisson equation
% 
%       -\Delta u - k^2u = f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition: (1)    g_R*u + d*grad(u)*n=g_N on \Gamma _R
%                                 (2)    -i*u + grad(u)*n = 0;   on \Gamma _R
%   The mesh is given by node and elem and the boundary edge is given by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, g_R, or d.
%   For general elliptic equations with convection and reaction
%   coefficients, see ellipticpde.
%
%   Example
%
%
%   See also Poisson3, squarePoisson, Lshape, crack, mg
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
N = size(node,1); NT = size(elem,1);
dim = size(node,2);
Ndof = N;

tic;  % record assembling time
%% Diffusion coefficient
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isempty(pde.d) && isnumeric(pde.d)
   K = pde.d;                                 % d is an array
end
if ~isempty(pde.d) && ~isnumeric(pde.d)       % d is a function   
    [lambda,weight] = quadpts(option.dquadorder);
    nQuad = size(lambda,1);
    K = zeros(NT,dim);
    for p = 1:nQuad
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
        temp = pde.d(pxy);
        K(:,1) = K(:,1) + weight(p)*temp(1:end/2);
        K(:,2) = K(:,2) + weight(p)*temp(1+end/2:end);
   end
end

%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
% Laplacian
Delta    = sparse(Ndof,Ndof);
M        = sparse(Ndof,Ndof);
for i = 1:3
    for j = i:3
        % $A_{ij}|_{\tau} = \int_{\tau}K\nabla \phi_i\cdot \nabla \phi_j dxdy$ 
        if ~isempty(pde.d)
            Aij = (Dphi(:,1,i).*K(:,1).*Dphi(:,1,j) + ...
                   Dphi(:,2,i).*K(:,2).*Dphi(:,2,j)).*area;
        else
            Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        end
        if (j==i)
            Delta = Delta + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
        else
            Delta = Delta + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                                   [Aij; Aij],Ndof,Ndof);        
        end 
        
        if ~isempty(pde.k2)
            fquadorder = 3;
             [lambda,weight] = quadpts(fquadorder);
                         phi = lambda;  % linear bases
	                   nQuad = size(lambda,1);
                         mij = 0;
            for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		mks1s2 = pde.k2(pxy);
            mij = mij + weight(p)*phi(p,i)*phi(p,j)*mks1s2;
            end
            mij = mij.*area;
        
        end
        if (j==i)
            M = M + sparse(elem(:,i),elem(:,j),mij,N,N);
        else    
            M = M + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [mij; mij],N,N);           
        end  
    end
end
% M;
% pause
clear K Aij
% mass matrix
%M0 = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
% if ~isempty(pde.k2) && isnumeric(pde.k2)
%    k2 = pde.k2;                           % k is an array
% end
% if ~isempty(pde.k2) && ~isnumeric(pde.k2)       % k is a function   
%    k2 = pde.k2(node);
% end

A = Delta - M;

%% The rest is the same as Poisson

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(area,1,3);
    b = accumarray(elem(:),bt(:),[Ndof 1]);
end
clear pxy bt

%% Set up boundary conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,b,u,freeNode] = getbd(b);




%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeNode), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        itStep = 1;
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'bicgstab'
                tic;
        [u(freeNode),~,~,iter] =  bicgstab(AD(freeNode,freeNode),b(freeNode),1.0e-6,N,Delta(freeNode,freeNode));
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'gmres'
                tic;
        %[u(freeNode),~,~,iter] =  gmres(AD(freeNode,freeNode),b(freeNode),[],1.0e-6,N);
        [u(freeNode),~,~,iter] =  gmres(AD(freeNode,freeNode),b(freeNode),[],1.0e-6,N,Delta(freeNode,freeNode));
        itStep = max(iter);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',itStep,'err',residual,'flag',2,'stopErr',residual);    
    case 'mg'
        option.x0 = u;
        option.solver = 'CG';
        [u,info] = mg(AD,b,elem,option);
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option);                 
end

%% Output information
eqn = struct('A',AD,'b',b,'freeNode',freeNode,'Delta',Delta);
info.assembleTime = assembleTime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeNode] = getbd(b)
    %% Set up of boundary conditions.
    %
    % 1) Modify the matrix for Dirichlet boundary nodes, which are not degree
    % of freedom. Values at these nodes are evaluatation of pde.g_D. The
    % original stiffness matrix A is turn into the matrix AD by enforcing
    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0.
    %
    % 2) Modify the right hand side b. The Neumann boundary integral is added
    % to b. For Dirichlet boundary ndoes, b(fixedNode) is the evaluation of
    % pde.g_D.
    %
    % Special attentation should be given for the pure Neumann boundary
    % condition. To enforce the compatible condition, the vector b should have
    % mean value zero. To avoid a singular matrix, the 1st node is chosen as
    % fixedNode. 
    %
    % The order of assigning Neumann and Dirichlet boundary condition is
    % important to get the right setting at the intersection nodes of Dirichlet
    % and Neumann boundary edges.
    %
    % Reference: Long Chen. Finite Element Methods and its Programming. Lecture
    % Notes.

    u = zeros(Ndof,1); 
    %% Initial check
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    %% Part 1: Modify the matrix for Dirichlet and Robin condition
    % Robin boundary condition
    Robin = [];
    isRobin = (bdFlag(:) == 3);
    if any(isRobin)
        allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
        Robin = allEdge(isRobin,:);
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        ve = node(Robin(:,1),:) - node(Robin(:,2),:);
        edgeLength = sqrt(sum(ve.^2,2)); 
        mid = (node(Robin(:,1),:) + node(Robin(:,2),:))/2;
        % int g_R phi_iphi_j ds
        ii = [Robin(:,1),Robin(:,1),Robin(:,2),Robin(:,2)];
        jj = [Robin(:,1),Robin(:,2),Robin(:,1),Robin(:,2)];
        temp = pde.g_R(mid).*edgeLength;
        ss = [1/3*temp, 1/6*temp, 1/6*temp, 1/3*temp];
        A = A + sparse(ii,jj,ss,Ndof,Ndof);
    end
    
    % Find Dirichlet boundary nodes: fixedNode
    fixedNode = []; freeNode = [];
    if ~isempty(bdFlag) % find boundary edges and boundary nodes
        [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
        freeNode = find(~isBdNode);
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given
        [fixedNode,bdEdge,isBdNode] = findboundary(elem);
        freeNode = find(~isBdNode);
    end
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0.
    if ~isempty(fixedNode)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedNode) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end
    
    %% Part 2: Find boundary edges and modify the right hand side b
    % Find boundary edges: Neumann
    Neumann = []; 
    if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
        Neumann = bdEdge;        
    end
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        [tempvar,Neumann] = findboundary(elem); %#ok<ASGLU>
    end

    % Neumann boundary condition
    if  isnumeric(pde.g_N) && all(pde.g_N == 0)
        pde.g_N = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.g_N)
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        phigN = lambdagN;                 % linear bases
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),2);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNp = pde.g_N(ppxy);
            for igN = 1:2
                ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
            end
        end
        ge = ge.*repmat(el,1,2);
        b = b + accumarray(Neumann(:), ge(:),[Ndof,1]); 
    end
    % The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann edges and no modification of
    % A,u,b is needed.

    % Dirichlet boundary condition
    if isnumeric(pde.g_D) && all(pde.g_D == 0)   % zero g_D
        pde.g_D = [];
    end
    if ~isempty(fixedNode) && ~isempty(pde.g_D)
        if isnumeric(pde.g_D)  % pde.g_D could be a numerical array 
            u(fixedNode) = pde.g_D(fixedNode); 
        else % pde.g_D is a function handle
            u(fixedNode) = pde.g_D(node(fixedNode,:));
        end
        b = b - A*u;
        b(fixedNode) = u(fixedNode);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.
    end % end of getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of Poisson