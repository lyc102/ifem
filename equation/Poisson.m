function [soln,eqn,info] = Poisson(node,elem,bdFlag,pde,option)
%% POISSON Poisson equation: P1 linear element.
%
%   u = Poisson(node,elem,bdFlag,pde) produces the linear finite element
%   approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
% 
%   The mesh is given by node and elem and the boundary edge is by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, g_R, or d.
%   For general elliptic equations with convection and reaction
%   coefficients, see ellipticpde.
%   
%   soln = Poisson(node,elem,bdFlag,pde,option) specifies the options.
%
%   In the output, soln structure contains
%     - soln.u: solution u
%     - soln.Du: gradient of u
%
%   In the input, option structures contains
%    - option.dquadorder: quadrature order for diffusion coefficients
%    - option.fquadorder: quadrature order for computing right hand side f
%    - option.solver
%      'direct': the built in direct solver \ (mldivide)
%      'mg':     multigrid-type solvers mg is used.
%      'amg':    algebraic multigrid method is used.
%      'none':   only assemble the matrix equation but not solve. 
%   The default setting is to use the direct solver for small size problems
%   and multigrid solvers for large size problems. For more options on the
%   multigrid solver mg, type help mg.
%
%   The function Poisson assembes the matrix equation AD*u = b and solves
%   it by the direct solver (small size <= 2e3) or the multigrid solver
%   (large size > 2e3). The Dirichlet boundary condition is built into the
%   matrix AD and the Neumann boundary condition is build into b.
%
%   The diffusion coefficient d is a scalar function or a column array with
%   the same length as the elem array. 
%
%   When only one type of boundary condition is imposed, the input argument
%   bdFlag can be skipped. The boundary condition is implicitly given in
%   the pde structure by specifying g_D or g_N only. See examples below.
%
%   [soln,eqn] = Poisson(node,elem,bdFlag,pde) returns also the equation
%   structure eqn, which includes: 
%     - eqn.AD:  modified stiffness matrix AD;
%     - eqn.b:   the right hand side. 
%     - eqn.Lap: non-modified stiffness matrix
%
%   The solution u = AD\b. The matrix eqn.Lap can be used to evulate the
%   bilinear form a(u,v) = u'*eqn.Lap*v, especially the enery norm of a finite
%   element function u is given by by sqrt(u'*eqn.Lap*u). 
%
%   [soln,eqn,info] = Poisson(node,elem,bdFlag,pde) returns also the
%   information on the assembeling and solver, which includes:
%     - info.assembleTime: time to assemble the matrix equation
%     - info.solverTime:   time to solve the matrix equation
%     - info.itStep:       number of iteration steps for the mg solver
%     - info.error:        l2 norm of the residual b - A*u
%     - info.flag:         flag for the mg solver.
%       flag = 0: converge within max iteration 
%       flag = 1: iterated maxIt times but did not converge
%       flag = 2: direct solver
%       flag = 3: no solve
%
%   Example
%     squarePoisson; Poissonfemrate;
%
%   Example
%     clear all
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:4
%       [node,elem] = uniformrefine(node,elem);
%     end
%     % Homogenous Dirichlet boundary condition
%     pde.f = inline('ones(size(p,1),1)','p');
%     pde.g_D = 0;
%     u = Poisson(node,elem,[],pde);
%     figure(1); 
%     showresult(node,elem,u);
%
%   See also Poisson3, femPoisson, squarePoisson, Poissonfemrate
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Preprocess
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
% important constants
N = size(node,1); 
NT = size(elem,1);
Ndof = N;

%% Diffusion coefficient
time = cputime;  % record assembling time
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isempty(pde.d) && isnumeric(pde.d)
   K = pde.d;                                 % d is an array
end
if ~isempty(pde.d) && ~isnumeric(pde.d)       % d is a function   
    [lambda,weight] = quadpts(option.dquadorder);
    nQuad = size(lambda,1);
    K = zeros(NT,1);
    for p = 1:nQuad
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
        K = K + weight(p)*pde.d(pxy);      
   end
end

%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);
for i = 1:3
    for j = i:3
        % $A_{ij}|_{\tau} = \int_{\tau}K\nabla \phi_i\cdot \nabla \phi_j dxdy$ 
        Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if ~isempty(pde.d)
            Aij = K.*Aij;
        end
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                           [Aij; Aij],Ndof,Ndof);        
        end        
    end
end
clear K Aij

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if isreal(pde.f) % f is a real number or vector and not a function
   switch length(pde.f)
       case NT  % f is piecewise constant
         bt = pde.f.*area/3;
         b = accumarray(elem(:),[bt; bt; bt],[Ndof 1]);
       case N   % f is piecewise linear
         bt = zeros(NT,3);
         bt(:,1) = area.*(2*pde.f(elem(:,1)) + pde.f(elem(:,2)) + pde.f(elem(:,3)))/12;
         bt(:,2) = area.*(2*pde.f(elem(:,2)) + pde.f(elem(:,3)) + pde.f(elem(:,1)))/12;
         bt(:,3) = area.*(2*pde.f(elem(:,3)) + pde.f(elem(:,1)) + pde.f(elem(:,2)))/12;
         b = accumarray(elem(:),bt(:),[Ndof 1]);
       case 1   % f is a scalar e.g. f = 1
         bt = pde.f*area/3;
         b = accumarray(elem(:),[bt; bt; bt],[Ndof 1]);
   end
end
if ~isempty(pde.f) && ~isreal(pde.f)  % f is a function 
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
[AD,b,u,freeNode,isPureNeumann] = getbd(b);

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeNode), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')  || isfield(option,'mgoption')   % no option.solver
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
        t = cputime;
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        if ~isfield(option,'mgoption')   % no option.mgoption
            option.mgoption.x0 = u;
            option.mgoption.solver = 'CG';
        end
        [u,info] = mg(AD,b,elem,option.mgoption);
    case 'amg'
        if ~isfield(option,'amgoption')  % no option.amgoption
            option.amgoption.x0 = u;
            option.amgoption.solver = 'CG';
        end
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option.amgoption);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    patchArea = accumarray(elem(:),[area;area;area]/3, [N 1]); 
    uc = sum(u.*patchArea)/sum(area);
    u = u - uc;   % int u = 0
end

%% Compute Du
dudx =  u(elem(:,1)).*Dphi(:,1,1) + u(elem(:,2)).*Dphi(:,1,2) ...
      + u(elem(:,3)).*Dphi(:,1,3);
dudy =  u(elem(:,1)).*Dphi(:,2,1) + u(elem(:,2)).*Dphi(:,2,2) ...
      + u(elem(:,3)).*Dphi(:,2,3);         
Du = [dudx, dudy];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'freeNode',freeNode,'Lap',A);
    info.assembleTime = assembleTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeNode,isPureNeumann] = getbd(b)
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
        % use Simplson rule to compute int g_R phi_iphi_j ds
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
    isPureNeumann = false;
    if isempty(fixedNode) && isempty(Robin) % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedNode = 1;
        freeNode = 2:Ndof;    % eliminate the kernel by enforcing u(1) = 0;
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
    if ~isPureNeumann && ~isempty(fixedNode) && ~isempty(pde.g_D)
        if isnumeric(pde.g_D)  % pde.g_D could be a numerical array 
            u(fixedNode) = pde.g_D(fixedNode); 
        else % pde.g_D is a function handle
            u(fixedNode) = pde.g_D(node(fixedNode,:));
        end
        b = b - A*u;
    end
    if ~isPureNeumann % non-empty Dirichlet boundary condition
        b(fixedNode) = u(fixedNode);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
        b(1) = 0;
    end
    end % end of getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of Poisson
