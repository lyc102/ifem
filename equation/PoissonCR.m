function [soln,eqn,info] = PoissonCR(node,elem,bdFlag,pde,option)
%% POISSONCR Poisson equation: Crouzeix-Raviart element.
%
%   u = POISSONCR(node,elem,bdFlag,pde) produces the Crouzeix and Raviart
%   linear nonconforming finite element approximation of the Poisson
%   equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
% [soln,eqn,info] = POISSONCR(node,elem,pde,bdFlag,option)
%
% The usage is the same as <a href="matlab:help Poisson">Poisson</a>. CR element on a triangle is
% summarized in <a href="matlab:ifem PoissonfemrateCR">PoissonfemrateCR</a> for detail.
% 
%   Example
%     clear all
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:3
%       [node,elem] = uniformrefine(node,elem);
%     end
%     % Homogenous Dirichlet boundary condition
%     pde.f = inline('ones(size(p,1),1)','p');
%     pde.g_D = 0;
%     u = PoissonCR(node,elem,[],pde);
%     figure(1); 
%     showresult(node,elem,u);
%
% Example
%
%    exampleCR
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Preprocess
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
% important constants
NT = size(elem,1);

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
[Dphi,area] = gradbasis(node,elem);
NE = size(edge,1); 
Ndof = NE;

%% Assemble stiffness matrix
A = sparse(Ndof,Ndof);
for i = 1:3
    for j = i:3
        % local to global index map
        ii = double(elem2edge(:,i));
        jj = double(elem2edge(:,j));
        % local stiffness matrix
        Aij = 4*dot(Dphi(:,:,i),Dphi(:,:,j),2).*area;
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
clear K Aij

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
    phi = 1-2*lambda;  % bases for CR element
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
    b = accumarray(elem2edge(:),bt(:),[Ndof 1]);
end
clear pxy bt

%% Set up boundary conditions
[AD,b,u,freeEdge,isPureNeumann] = getbdCR(b);

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeEdge), return; end
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
    u(freeEdge) = AD(freeEdge,freeEdge)\b(freeEdge);         
    residual = norm(b - AD*u);
    info = struct('solverTime',toc,'itStep',0,'error',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'error',[],'flag',3,'stopErr',[]);
    case 'mg'
        option.x0 = u;
        option.solver = 'CG';
        [u,info] = mg(AD,b,elem,option,edge);
    case 'amg'
        option.solver = 'CG';
        [u(freeEdge),info] = amg(AD(freeEdge,freeEdge),b(freeEdge),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    patchArea = accumarray(elem2edge(:),[area;area;area]/3, [NE 1]);     
    uc = sum(u.*patchArea)/sum(area);
    u = u - uc;   % normalization for pure Neumann problem
end

%% Compute Du
dudx =  u(elem2edge(:,1)).*Dphi(:,1,1) + u(elem2edge(:,2)).*Dphi(:,1,2) ...
      + u(elem2edge(:,3)).*Dphi(:,1,3);
dudy =  u(elem2edge(:,1)).*Dphi(:,2,1) + u(elem2edge(:,2)).*Dphi(:,2,2) ...
      + u(elem2edge(:,3)).*Dphi(:,2,3);         
Du = -2*[dudx, dudy];

%% Output information
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'edge',edge,'freeEdge',freeEdge);
    info.assembleTime = assembleTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeEdge,isPureNeumann]= getbdCR(b)
    %% GETBDCR Boundary conditions for Poisson equation: Crouzeix-Raviart element.
    
    u = zeros(Ndof,1);
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
        ii = find(isRobin);
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
        s = accumarray(elem2edge(:), 1, [Ndof 1]);
        fixedEdge = find(s == 1);
        freeEdge = find(s == 2);
    end
    isPureNeumann = false;    
    if isempty(fixedEdge) && isempty(Robin)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedEdge = 1;
        freeEdge = (2:Ndof)';    % eliminate the kernel by enforcing u(1) = 0;
    end
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedEdge,fixedEdge)=I, AD(fixedEdge,freeEdge)=0, AD(freeEdge,fixedEdge)=0.
    if ~isempty(fixedEdge)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedEdge) = 1;
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
        b(isNeumann) = b(isNeumann) + ge;
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
            u(fixedEdge) = pde.g_D(mid);
        end
        b = b - A*u;
    end
    if ~isPureNeumann % non-empty Dirichlet boundary condition
        b(fixedEdge) = u(fixedEdge);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
        b(1) = 0;
    end
    end % end of getbdCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of PoissonCR