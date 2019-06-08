function [soln,eqn,info] = PoissonP2(node,elem,bdFlag,pde,option)
%% POISSONP2 Poisson equation: P2 quadratic element.
%
% u = PoissonP2(node,elem,bdFlag,pde,option) produces the quadratic
%   finite element approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
% [soln,eqn,info] = PoissonP2(node,elem,pde,bdFlag,option)
%
% The usage is the same as <a href="matlab:help Poisson">Poisson</a>. Quadratic element on a triangle is
% summarized in <a href="matlab:ifem PoissonP2femrate">PoissonP2femrate</a> for detail.
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
%     u = PoissonP2(node,elem,[],pde);
%     figure(1); 
%     showresult(node,elem,u);
%
% See also Poisson, Poisson3, Poisson3P2 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Preprocess
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end

%% Construct data structure
time = cputime;  % record assembling time
[elem2dof,edge,bdDof] = dofP2(elem);
% important constants
N = size(node,1);  NT = size(elem,1); NE = size(edge,1);
Ndof = N + NE;
elem2edge = elem2dof(:,4:6) - N;

%% Compute geometric quantities and gradient of local basis
[Dlambda,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
% generate sparse pattern
ii = zeros(21*NT,1); 
jj = zeros(21*NT,1); 
index = 0;
for i = 1:6
    for j = i:6
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));  
        index = index + NT;
    end
end
% quadrature points
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'quadorder')
    % diffusion is piecewise constant
    option.quadorder = 2;        % default order
    if ~isempty(pde.d) && isnumeric(pde.d) % numerical diffusion
        option.quadorder = 3;    % exact for linear diffusion coefficient
    end
end
[lambda, w] = quadpts(option.quadorder);
nQuad = size(lambda,1);
% compute non-zeros
sA = zeros(21*NT,nQuad);
for p = 1:nQuad
    % Dphi at quadrature points
    Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
    Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
    Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
    Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
    Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
    index = 0;
    for i = 1:6
        for j = i:6
            Aij = 0;
            if isempty(pde.d) || isnumeric(pde.d)
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            else
                pxy = lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:);
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.d(pxy);
            end
            if ~isempty(pde.d) && isnumeric(pde.d) % d is piecewise constant
                Aij = pde.d.*Aij;
            end
            Aij = Aij.*area;
            sA(index+1:index+NT,p) = Aij;
            index = index + NT;
        end
    end
end
sA = sum(sA,2);
% assemble the matrix
diagIdx = (ii == jj);   
upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj sA

%% Assemble right hand side by high order quadrature rule
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 4;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    phi(:,6) = 4*lambda(:,1).*lambda(:,2);
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phi(:,4) = 4*lambda(:,2).*lambda(:,3);
    phi(:,5) = 4*lambda(:,3).*lambda(:,1);
    bt = zeros(NT,6);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        if isfield(pde,'f') && isnumeric(pde.f)
            fp = pde.f;        % piecewise constant       
        else
            fp = pde.f(pxy);   % function handle
        end
        for j = 1:6
            bt(:,j) = bt(:,j) + w(p)*phi(p,j)*fp;
        end
    end
    bt = bt.*repmat(area,1,6);
    b = accumarray(elem2dof(:),bt(:),[Ndof 1]); 
end

%% Boundary Conditions
[AD,b,u,freeDof,isPureNeumann] = getbdP2(b);

%% Record assembeling time
assembleTime = cputime - time;
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
    else            % Multigrid-type  solver for large size systems
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
        option.x0 = u;
        option.solver = 'CG';
        option.tol = 1e-9;
        [u,info] = mg(AD,b,elem,option,edge);
    case 'amg'
        option.solver = 'CG';
        option.tol = 1e-9;
        [u(freeDof),info] = amg(AD(freeDof,freeDof),b(freeDof),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
      patchArea = accumarray(elem2edge(:),[area;area;area]/3, [NE 1]);     
             uc = sum(u(N+1:end).*patchArea)/sum(area);
              u = u - uc;   % normalization for pure Neumann problem     
% Here we use 3 middle points rule to compute the integral which is exact
% for quadratic elements. We can also use the following formulae     
%      intuh = sum(1/3*u(elem2dof(:,4:6)),2).*area; %Compute the integrable of uh
%         uc = sum(intuh)/sum(area);
%          u = u -uc;
end

%% Compute Du
Du = [];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'edge',edge,'freeDof',freeDof,'Lap',A);
    info.assembleTime = assembleTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdP2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeDof,isPureNeumann] = getbdP2(b)
    %% Boundary conditions for Poisson equation: P2 quadratic FEM.
    %
    % The set up of boundary condition consists of two parts: 
    %
    % 1) Modify the matrix for Dirichlet boundary nodes, which are not degree
    % of freedom. Values at these nodes are evaluatation of pde.g_D. The
    % original stiffness matrix A is turn into the matrix AD by enforcing
    % AD(fixedDof,fixedDof)=I, AD(fixedDof,freeDof)=0, AD(freeDof,fixedDof)=0.
    %
    % 2) Modify the right hand side b. The Neumann boundary integral is added
    % to b. For Dirichlet boundary ndoes, b(fixedDof) is the evaluation of
    % pde.g_D.
    %
    % Special attentation should be given for the pure Neumann boundary
    % condition. To enforce the compatible condition, the vector b should have
    % mean value zero. To avoid a singular matrix, the 1st node is chosen as
    % fixedDof. 
    %
    % The order of assigning Neumann and Dirichlet boundary condition is
    % important to get the right setting at the intersection nodes of Dirichlet
    % and Neumann boundary edges.

    u = zeros(Ndof,1);
    %% Set up boundary and basic parameter
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
        if ~isfield(option,'gRquadorder')
            option.gRquadorder = 5;   % default order exact for linear gR
        end
        [lambdagR,weightgR] = quadpts1(option.gRquadorder);
        nQuadgR = size(lambdagR,1);
        % quadratic bases (1---3---2)
        bdphi = zeros(nQuadgR,3);
        bdphi(:,1) = (2*lambdagR(:,1)-1).*lambdagR(:,1);
        bdphi(:,2) = (2*lambdagR(:,2)-1).*lambdagR(:,2);
        bdphi(:,3) = 4*lambdagR(:,1).*lambdagR(:,2);
        % length of edge
        el = sqrt(sum((node(Robin(:,1),:) - node(Robin(:,2),:)).^2,2));
        NR = size(Robin,1);
        ss = zeros(NR,3,3);
        for pp = 1:nQuadgR
            ppxy = lambdagR(pp,1)*node(Robin(:,1),:) ...
                 + lambdagR(pp,2)*node(Robin(:,2),:);
            gRp = pde.g_R(ppxy);
            for iR = 1:3
                for jR = iR:3   % only compute half of the off-diagonal part
                    ss(:,iR,jR) = ss(:,iR,jR) + ...
                    weightgR(pp)*gRp*bdphi(pp,iR).*bdphi(pp,jR);
                end
            end
        end
        ss(:) = ss(:).*repmat(el,9,1);
        index = 0;
        Robin(:,3) = find(isRobin)+N;  % the third one maps to corresponding dof
        for iR = 1:3
            for jR = 1:3
                iiR(index+1:index+NR) = double(Robin(:,iR)); 
                jjR(index+1:index+NR) = double(Robin(:,jR)); 
                if jR>=iR
                    ssR(index+1:index+NR) = ss(:,iR,jR);
                else
                    ssR(index+1:index+NR) = ss(:,jR,iR);
                end
                index = index + NR;
            end
        end
        A = A + sparse(iiR,jjR,ssR,Ndof,Ndof);
    end

    % Find Dirichlet boundary dof: fixedDof
    fixedDof = []; 
    freeDof = [];
    isFixedDof = false(Ndof,1); 
    if ~isempty(bdFlag)     
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
        isFixedDof(edge(isDirichlet,:)) = true;
        isFixedDof(N + find(isDirichlet')) = true;
        fixedDof = find(isFixedDof);
        freeDof = find(~isFixedDof);    
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        fixedDof = bdDof;
        isFixedDof(fixedDof) = true;
        freeDof = find(~isFixedDof);    
    end
    isPureNeumann = false;        
    if isempty(fixedDof) && isempty(Robin)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedDof = 1;
        freeDof = 2:Ndof;    % eliminate the kernel by enforcing u(1) = 0;
    end
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % |AD(fixedDof,fixedDof)=I, AD(fixedDof,freeDof)=0,
    % AD(freeDof,fixedDof)=0|.
    if ~isempty(fixedDof)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedDof) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end

    %% Part 2: Find boundary edges and modify the load b
    % Find boundary edges: Neumann and Robin
    Neumann = [];
    if ~isempty(bdFlag)     
        idxN = (bdFlag(:) == 2);      % all Neumann edges in bdFlag        
        Neumannidx = elem2edge(idxN | idxR); % index of Neumann and Robin edges
        % since boundary integral is also needed for Robin edges
        Neumann   = edge(Neumannidx,:);
    end
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        Neumannidx = find(bdDof>N);
        Neumann = edge(Neumannidx,:);
    end

    % Neumann boundary condition
    if ~isempty(pde.g_N) && ~isempty(Neumann) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 4;   % default order
        end
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        nQuadgN = size(lambdagN,1);
        % quadratic bases on an edge (1---3---2)
        bdphi = zeros(nQuadgN,3);        
        bdphi(:,1) = (2*lambdagN(:,1)-1).*lambdagN(:,1);
        bdphi(:,2) = (2*lambdagN(:,2)-1).*lambdagN(:,2);
        bdphi(:,3) = 4*lambdagN(:,1).*lambdagN(:,2);
        % length of edge
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        ge = zeros(size(Neumann,1),3);
        for pp = 1:nQuadgN
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNp = pde.g_N(ppxy);
            ge(:,1) = ge(:,1) + weightgN(pp)*gNp*bdphi(pp,1);
            ge(:,2) = ge(:,2) + weightgN(pp)*gNp*bdphi(pp,2);
            ge(:,3) = ge(:,3) + weightgN(pp)*gNp*bdphi(pp,3); % interior bubble
        end
        % update RHS
        ge = ge.*repmat(el,1,3);        
        b(1:N) = b(1:N) + accumarray(Neumann(:), [ge(:,1); ge(:,2)],[N,1]);
        b(N+Neumannidx) = b(N+Neumannidx) + ge(:,3);
    end

    % Dirichlet boundary conditions
    if ~isPureNeumann && ~isempty(fixedDof) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
        idx = (fixedDof > N);                         % index of edge nodes
        u(fixedDof(~idx)) = pde.g_D(node(fixedDof(~idx),:)); % bd value at vertex dofs
        bdEdgeIdx = fixedDof(idx) - N;
        bdEdgeMid = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2;
        u(fixedDof(idx)) = pde.g_D(bdEdgeMid);
        b = b - A*u;
    end
    if ~isPureNeumann % non-empty Dirichlet boundary condition
        b(fixedDof) = u(fixedDof);
    end

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
        b(1) = 0;        
    end
    end % end of getbdP2
end % end of function PoissonP2