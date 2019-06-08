function [soln,eqn,info] = Poisson3WG(node,elem,bdFlag,pde,option,varargin)
%% POISSON3WG Poisson equation: P1 linear element in 3-D.
%
%   u = POISSON3WG(node,elem,pde,bdFlag) produces the linear finite element
%   approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
% 
%   The mesh is given by node and elem and the boundary face is given by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, g_R, or d.
%   For general elliptic equations with convection and reaction
%   coefficients, see ellipticpde.
%   
%   The function Poisson assembes the matrix equation AD*u = b and solves
%   it by the direct solver (small size <= 2e3) or the multigrid solver
%   (large size > 2e3). The Dirichlet boundary condition is built into the
%   matrix AD and the Neumann boundary condition is build into b.
%
%   The diffusion coefficient d is a scalar function or a column array with
%   the same length as the elem array. 
% 
%   u = Poisson3(node,elem,pde,bdFlag,HB,option) specifies the options.
%     - option.solver == 'direct': the built in direct solver \ (mldivide)
%     - option.solver == 'mg':     multigrid-type solvers mg is used.
%     - option.solver == 'notsolve': the solution u = u_D. 
%   The default setting is to use the direct solver for small size problems
%   and multigrid solvers for large size problems. For more options on the
%   multigrid solver mg, type help mg.
%
%   When only one type of boundary condition is imposed, the input argument
%   bdFlag can be skipped. The boundary condition is implicitly given in
%   the pde structure by specifying g_D or g_N only. See examples below.
%
%   [u,A] = Poisson3(node,elem,pde,bdFlag,HB) returns also the
%   non-modified stiffness matrix A, which is semi-definite. The kernel of
%   A consists of constant vectors. The matrix A can be used to evulate the
%   bilinear form A(u,v) = u'*A*v, especially the enery norm of a finite
%   element function u by sqrt(u'*A*u).
%
%   [u,A,eqn] = Poisson3(node,elem,pde,bdFlag,HB) returns also the equation
%   structure eqn, which includes: 
%     - eqn.AD:  modified stiffness matrix AD;
%     - eqn.b:   the right hand side. 
%   The solution u = AD\b. The output eqn can be used to test other solvers.
%
%   [u,A,eqn,info] = Poisson3(node,elem,pde,bdFlag,HB) returns also the
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
%     clear all
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:4
%       [node,elem] = uniformbisect(node,elem);
%     end
%     % Homogenous Dirichlet boundary condition
%     pde.f = inline('ones(size(p,1),1)','p');
%     pde.g_D = 0;
%     u = Poisson(node,elem,pde);
%     figure(1); 
%     showresult(node,elem,u);
%     % Non-homogenous Dirichlet boundary condition
%     pde.f = inline('-4*ones(size(p,1),1)','p');
%     pde.g_D = inline('p(:,1).^2 + p(:,2).^2','p');
%     u = Poisson(node,elem,pde);
%     figure(2); 
%     showresult(node,elem,u);
%     % Homogenous Neumann boundary condition
%     clear pde
%     pde.f = inline('pi^2*cos(pi*p(:,1)).*cos(pi*p(:,2))','p');
%     u = Poisson(node,elem,pde);
%     figure(3);
%     showresult(node,elem,u);
%
%   Example
%     cubePoisson;
%
%   See also Poisson3, squarePoisson, Lshape, crack, mg
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.


if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
NT = size(elem,1);

%% Diffusion coefficient
t = cputime;  % record assembling time
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isempty(pde.d) && isnumeric(pde.d)
   K = pde.d;                   % d is an array
end
if ~isempty(pde.d) && ~isnumeric(pde.d)
    [lambda,weight] = quadpts3(option.dquadorder);
    nQuad = size(lambda,1);
    K = zeros(NT,1);
    for p = 1:nQuad
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:) ...
            + lambda(p,4)*node(elem(:,4),:);
        K = K + weight(p)*pde.d(pxy);           % d is a function   
   end
end

%% Compute geometric quantities and gradient of local basis
[elem2face,face] = dof3face(elem);
NF = size(face,1); 
Ndof = NT + NF;
elem2dof = NT + elem2face;

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);
% compute ct2 = 1/mean(||x-xc||^2)
center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
          node(elem(:,3),:) + node(elem(:,4),:))/4;
[lambda,weight] = quadpts3(2);
nQuad = size(lambda,1);
ct2 = zeros(NT,1);
for p = 1:nQuad
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
         + lambda(p,2)*node(elem(:,2),:) ...
         + lambda(p,3)*node(elem(:,3),:) ...
         + lambda(p,4)*node(elem(:,4),:);
    ct2 = ct2 + weight(p)*sum((pxyz-center).^2,2); 
end
ct2 = 1./ct2;
[Dphi,volume] = gradbasis3(node,elem);
clear center pxyz

% Mbb: face - face           
for i = 1:4
    for j = i:4
        % local to global index map
		ii = double(elem2dof(:,i));
		jj = double(elem2dof(:,j)); 
        % compuation of local stiffness matrix.        
        Aij = 9*dot(Dphi(:,:,i),Dphi(:,:,j),2).*volume + 9/16*ct2.*volume;
        if ~isempty(pde.d)
            Aij = K.*Aij; 
        end 
        if (j==i)
            A = A + sparse(ii,jj,Aij,Ndof,Ndof);
        else
            A = A + sparse([ii,jj],[jj,ii],[Aij; Aij],Ndof,Ndof);        
        end        
    end
end

% Mob: interior - face
Aij = -9/4*ct2.*volume;
if ~isempty(pde.d)
    Aij = K.*Aij;
end
Mob = sparse([(1:NT)', (1:NT)', (1:NT)', (1:NT)'], ...
             double(elem2dof(:)), [Aij, Aij, Aij, Aij], Ndof, Ndof);
A = A + Mob + Mob';

% Moo: diagonal of interor
Aij = 9*ct2.*volume;
if ~isempty(pde.d)
    Aij = K.*Aij;
end
A =  A + sparse(1:NT, 1:NT, Aij, Ndof, Ndof);

clear K Aij

%% Assemble right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    [lambda,weight] = quadpts3(option.fquadorder);
	nQuad = size(lambda,1);
    bt = zeros(NT,1);
    for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxyz);
        bt = bt + weight(p)*fp;
    end
    bt = bt.*volume;
    b(1:NT) = bt;
end
clear pxyz bt

%% Set up boundary conditions
[AD,b,u,freeDof,isPureNeumann] = getbdWG3(A,b);

%% Record assembling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeDof), return; end
% Set up solver
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
        u(freeDof) = AD(freeDof,freeDof)\b(freeDof);         
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'error',[],'flag',3);
    case 'mg'
        option.solver = 'CG';
        if nargin>=6
            HB = varargin{1};
        else
            HB = [];
        end        
        if isfield(option,'reducesystem') && (option.reducesystem == 0)
            % original system
            option.x0 = u;
            option.freeDof = freeDof;
            [u,info] = mg(AD,b,elem,option,face,HB);        
        else  % reduced system
            option.x0 = u(NT+1:end); 
            option.freeDof = freeDof(NT+1:end)-NT;
            % eleminate elementwise dof
            Aoinv = spdiags(1./diag(AD(1:NT,1:NT)),0,NT,NT);
            Aob = AD(1:NT,NT+1:end);
            Abo = Aob';
            Abb = AD(NT+1:end,NT+1:end);
            Abbm = Abb - Abo*Aoinv*Aob;
            bm = -Abo*Aoinv*b(1:NT) + b(NT+1:end);        
            [ub,info] = mg(Abbm,bm,elem,option,face,HB);
            u(1:NT) = Aoinv*(b(1:NT) - Aob*ub);
            u(NT+1:end) = ub;        
        end
    case 'amg'
        option.solver = 'CG';
        [u(freeDof),info] = amg(AD(freeDof,freeDof),b(freeDof),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    uc = sum(u(1:NT).*volume);
    u = u - uc;   % normalization for pure Neumann problem
end

%% Compute Du
dudx = u(elem2dof(:,1)).*Dphi(:,1,1)+u(elem2dof(:,2)).*Dphi(:,1,2) ...
     + u(elem2dof(:,3)).*Dphi(:,1,3)+u(elem2dof(:,4)).*Dphi(:,1,4);
dudy = u(elem2dof(:,1)).*Dphi(:,2,1)+u(elem2dof(:,2)).*Dphi(:,2,2) ...
     + u(elem2dof(:,3)).*Dphi(:,2,3)+u(elem2dof(:,4)).*Dphi(:,2,4);
dudz = u(elem2dof(:,1)).*Dphi(:,3,1)+u(elem2dof(:,2)).*Dphi(:,3,2) ...
     + u(elem2dof(:,3)).*Dphi(:,3,3)+u(elem2dof(:,4)).*Dphi(:,3,4);
Du = -3*[dudx, dudy, dudz]; 

%% Output information
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'face',face,'freeDof',freeDof,'Lap',A);
    info.assembleTime = assembleTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,b,u,freeDof,isPureNeumann]= getbdWG3(A,b)
    %% GETBDCR Boundary conditions for Poisson equation: Crouzeix-Raviart element.
    
    u = zeros(Ndof,1);
    %% Initial check
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    %% Part 1: Modify the matrix for Dirichlet and Robin condition
    % Robin boundary condition
    Robin = [];
    idxR = (bdFlag(:) == 3);      % index of Robin faces in bdFlag
    if any(idxR)    
        isRobin = false(Ndof,1);
        isRobin(elem2face(idxR)) = true;
        Robin = face(isRobin,:);  % Robin faces  
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        v12 = node(Robin(:,2),:)-node(Robin(:,1),:);
        v13 = node(Robin(:,3),:)-node(Robin(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        center = (node(Robin(:,1),:) + node(Robin(:,2),:)+node(Robin(:,3),:))/3;
        ii = NT + find(isRobin);
        ss = pde.g_R(center).*area; % exact for linear g_R
        A = A + sparse(ii,ii,ss,Ndof,Ndof);
    end
    
    % Dirichlet boundary nodes: fixedFace
    fixedFace = []; freeFace= [];
    if ~isempty(bdFlag) % find boundary faces
        idxD = (bdFlag(:)==1); % all Dirichlet faces in bdFlag
        isFixedFace = false(NF,1);
        isFixedFace(elem2face(idxD)) = true;
        fixedFace = find(isFixedFace);
        freeFace = find(~isFixedFace);
    end    
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given in the input
        s = accumarray(elem2face(:), 1, [NF 1]);
        fixedFace = find(s == 1);
        freeFace = find(s == 2);
    end     
    isPureNeumann = false;    
    if isempty(fixedFace) && isempty(Robin)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedFace = 1;
        freeFace = (2:NF)';    % eliminate the kernel by enforcing u(1) = 0;
    end    
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedFace,fixedFace)=I, AD(fixedFace,freeFace)=0, AD(freeFace,fixedFace)=0.
    if ~isempty(fixedFace)
        bdidx = zeros(Ndof,1); 
        bdidx(NT + fixedFace) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end
    
    %% Part 2: Find boundary faces and modify the right hand side b
    % Find boundary faces: Neumann
    Neumann = [];
    if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
        idxN = (bdFlag(:) == 2);      % all Neumann faces in bdFlag
        isNeumann = elem2face(idxN | idxR); % index of Neumann faces
        % since boundary integral is also needed for Robin edges
        Neumann = face(isNeumann,:);      % Neumann faces        
    end    
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        s = accumarray(elem2face(:), 1, [NF 1]);
        Neumann = face(s == 1,:);
    end
    
    % Neumann boundary condition
    if ~isempty(Neumann) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 3;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts(option.gNquadorder);
        nQuadgN = size(lambdagN,1);
        gf = zeros(size(Neumann,1),1);
        v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
        v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:)...
                 + lambdagN(pp,3)*node(Neumann(:,3),:);
            gNp = pde.g_N(ppxy);
            gf = gf+ weightgN(pp)*gNp;
        end
        gf = gf.*area;
        b(NT+isNeumann) = b(NT+isNeumann) + gf;
    end        
    % The case with non-empty Neumann faces but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann faces and no modification of
    % A,u,b is needed.
   
    % Dirichlet boundary condition
    if ~isPureNeumann && ~isempty(fixedFace) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
        if isnumeric(pde.g_D)
            u(fixedFace) = pde.g_D(fixedFace);
        else % pde.g_D is a function handle
            center = (node(face(fixedFace,1),:)+node(face(fixedFace,2),:)+node(face(fixedFace,3),:))/3;
            u(NT+fixedFace) = pde.g_D(center);
        end
        b = b - A*u;
        b(NT+fixedFace) = u(NT+fixedFace);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b); % compatilbe condition: sum(b) = 0
        b(1) = 0;
    end
    
    freeDof = [(1:NT)'; NT+freeFace];    
    end % end of getbd3WG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end