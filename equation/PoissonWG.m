function [soln,eqn,info] = PoissonWG(node,elem,bdFlag,pde,option)
%% POISSONWG Poisson equation: lowest order weak Galerkin element
%
%   u = POISSONWG(node,elem,bdFlag,pde) produces the linear finite element
%   approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
% 
% The usage is the same as <a href="matlab:help Poisson">Poisson</a>. Weak Galerkin method on a triangle is
% summarized in <a href="matlab:ifem PoissonWGfemrate">PoissonWGfemrate</a> for detail.
%
%   Example
%
%     squarePoissonWG;
%
%   See also Poisson, squarePoissonWG
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Preprocess
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
% important constants
NT = size(elem,1); 
N = size(node,1);

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

%% Construct data structure 
[elem2edge,edge] = dofedge(elem);
NE = size(edge,1); 
Ndof = NT + NE;
elem2dof = NT + elem2edge;

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);

% compute ct2 = 1/mean(||x-xc||^2)
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
mid1 = (node(elem(:,2),:) + node(elem(:,3),:))/2;
mid2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;
mid3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;
ct2 = 3./sum((mid1 - center).^2 + (mid2 - center).^2 + (mid3 - center).^2,2);
[Dphi,area] = gradbasis(node,elem);
clear center mid1 mid2 mid3

% Mbb: edge - edge           
for i = 1:3
    for j = i:3
        % local to global index map
        ii = double(elem2dof(:,i));
        jj = double(elem2dof(:,j));
        % local stiffness matrix
        Aij = 4*dot(Dphi(:,:,i),Dphi(:,:,j),2).*area + 4/9*ct2.*area;
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

% Mob: interior - edge
Aij = -4/3*ct2.*area;
if ~isempty(pde.d)
    Aij = K.*Aij;
end
Mob = sparse([(1:NT)', (1:NT)', (1:NT)'], ...
             double(elem2dof(:)), [Aij, Aij, Aij], Ndof, Ndof);
A = A + Mob + Mob';

% Moo: diagonal of interor
Aij = 4*ct2.*area;
if ~isempty(pde.d)
    Aij = K.*Aij;
end
A =  A + sparse(1:NT, 1:NT, Aij, Ndof, Ndof);

clear K Aij

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if isreal(pde.f) % f is a real number or vector and not a function
   switch length(pde.f)
       case NT  % f is piecewise constant
         b(1:NT) = pde.f.*area;
       case N   % f is piecewise linear
         b(1:NT) = (pde.f(elem(:,1)) + pde.f(elem(:,2)) + pde.f(elem(:,3)))/3.*area;
       case 1   % f is a scalar e.g. f = 1
         b(1:NT) = pde.f*area;
   end
end
if ~isempty(pde.f) && ~isreal(pde.f)  % f is a function
    [lambda,weight] = quadpts(option.fquadorder);
	nQuad = size(lambda,1);
    bt = zeros(NT,1);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
        bt = bt + weight(p)*fp;
    end
    bt = bt.*area;
    b(1:NT) = bt;
end
clear pxy bt

%% Set up boundary conditions
if nargin<=3, bdFlag = []; end
[AD,b,u,freeDof,isPureNeumann] = getbdWG(A,b);

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeDof), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if NE <= 1e3  % Direct solver for small size systems
        option.solver = 'direct';
    else          % MGCG  solver for large size systems
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
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        % eleminate elementwise dof
        option.solver = 'CG';
        if isfield(option,'reducesystem') && (option.reducesystem == 0)
            option.x0 = u;
            [u,info] = mg(AD,b,elem,option,edge);        
        else
            option.x0 = u(NT+1:end);
            Aoinv = spdiags(1./diag(AD(1:NT,1:NT)),0,NT,NT);
            Aob = AD(1:NT,NT+1:end);
            Abo = Aob';
            Abb = AD(NT+1:end,NT+1:end);
            Abbm = Abb - Abo*Aoinv*Aob;
            bm = -Abo*Aoinv*b(1:NT) + b(NT+1:end);        
            [ub,info] = mg(Abbm,bm,elem,option,edge);
            u(1:NT) = Aoinv*(b(1:NT) - Aob*ub);
            u(NT+1:end) = ub;        
        end
    case 'amg'
        option.solver = 'CG';
        [u(freeDof),info] = amg(AD(freeDof,freeDof),b(freeDof),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    uc = sum(u(1:NT).*area);
    u = u - uc;   % normalization for pure Neumann problem
end

%% Compute Du
dudx =  u(elem2dof(:,1)).*Dphi(:,1,1) + u(elem2dof(:,2)).*Dphi(:,1,2) ...
      + u(elem2dof(:,3)).*Dphi(:,1,3);
dudy =  u(elem2dof(:,1)).*Dphi(:,2,1) + u(elem2dof(:,2)).*Dphi(:,2,2) ...
      + u(elem2dof(:,3)).*Dphi(:,2,3);         
Du = -2*[dudx, dudy];

%% Output information
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'edge',edge,'freeDof',freeDof);
    info.assembleTime = assembleTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdWG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,b,u,freeDof,isPureNeumann]= getbdWG(A,b)
    %% GETBDCR Boundary conditions for Poisson equation: WG element.
    
    u =zeros(Ndof,1);
    %% Initial check
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    %% Part 1: Modify the matrix for Dirichlet and Robin condition
    % Robin boundary condition
    Robin = [];
    idxR = (bdFlag(:) == 3);      % index of Robin edges in bdFlag
    if any(idxR)    
        isRobin = false(NE,1);
        isRobin(elem2edge(idxR)) = true;
        Robin = edge(isRobin,:);  % Robin edges  
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        ve = node(Robin(:,1),:) - node(Robin(:,2),:);
        edgeLength = sqrt(sum(ve.^2,2)); 
        mid = (node(Robin(:,1),:) + node(Robin(:,2),:))/2;
        ii = NT + find(isRobin);  % for WG: edge dof is after elem dof
        ss = pde.g_R(mid).*edgeLength; % exact for linear g_R
        A = A + sparse(ii,ii,ss,Ndof,Ndof);
    end
    
    % Find Dirichlet boundary nodes: fixedEdge
    fixedEdge = []; freeEdge = [];
    if ~isempty(bdFlag)              % find boundary edges
        idxD = (bdFlag(:) == 1);     % all Dirichlet edges in bdFlag
        isFixedEdge = false(NE,1);
        isFixedEdge(elem2edge(idxD)) = true;  % index of fixed boundary edges
        fixedEdge = find(isFixedEdge);
        freeEdge = find(~isFixedEdge);
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given in the input
        s = accumarray(elem2edge(:), 1, [NE 1]);
        fixedEdge = find(s == 1);
        freeEdge = find(s == 2);
    end
    isPureNeumann = false;    
    if isempty(fixedEdge) && isempty(Robin)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedEdge = 1;
        freeEdge = (2:NE)';    % eliminate the kernel by enforcing u(1) = 0;
    end
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedEdge,fixedEdge)=I, AD(fixedEdge,freeEdge)=0, AD(freeEdge,fixedEdge)=0.
    if ~isempty(fixedEdge)
        bdidx = zeros(Ndof,1); 
        bdidx(NT + fixedEdge) = 1;
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
        idxN = (bdFlag(:) == 2);      % all Neumann edges in bdFlag
        isNeumann = elem2edge(idxN | idxR); % index of Neumann and Robin edges
        % since boundary integral is also needed for Robin edges
        Neumann = edge(isNeumann,:);      % Neumann edges        
    end
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        s = accumarray(elem2edge(:), 1, [NE 1]);
        Neumann = edge(s == 1,:);
    end

    % Neumann boundary condition
    if ~isempty(Neumann) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 3;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),1);
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNp = pde.g_N(ppxy);
            ge = ge+ weightgN(pp)*gNp;
        end
        ge = ge.*el;
        b(NT+isNeumann) = b(NT+isNeumann) + ge;
    end
    % The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann edges and no modification of
    % A,u,b is needed.

    % Dirichlet boundary condition
    if ~isPureNeumann && ~isempty(fixedEdge) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
        if isnumeric(pde.g_D)  % pde.g_D could be a numerical array 
            u(fixedEdge) = pde.g_D(fixedEdge); 
        else % pde.g_D is a function handle
            mid = (node(edge(fixedEdge,1),:) + node(edge(fixedEdge,2),:))/2;
            u(NT+fixedEdge) = pde.g_D(mid);
        end
        b = b - A*u;
        b(NT+fixedEdge) = u(NT+fixedEdge);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
        b(1) = 0;
    end
    
    freeDof = [(1:NT)'; NT+freeEdge];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of PoissonWG
