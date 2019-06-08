function [sigma,u,eqn,info] = HodgeLaplacianE(node,elem,bdFlag,pde,option)
%% HODEGELAPLACIANE Hodge Laplacian of edge element in two dimensions
%
% [sigma,u] = HodgeLaplacianE(node,elem,bdFlag,pde) computes the mixed
% finite element method approximation of the Hodge Laplacian equation
%
% Multigrid solvers are added by Jie Zhou.

if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

t = cputime;
%% Data structure
elemMG = elem;     % save elem and bdFlag for multigrid
bdFlagMG = bdFlag;
[elem,bdFlag] = sortelem(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dofedge(elem);
[Dlambda,area,elemSign] = gradbasis(node,elem);
locEdge = [2 3; 1 3; 1 2];
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
Nsigma = N; Nu = NE; Ndof = Nsigma + Nu;

%% Assemble matrix 
% Mass matrices
Mv = getmassmat(node,elem,area);
Me = getmassmatvec(elem2edge,area,Dlambda,'ND0');
invMt = spdiags(1./area,0,NT,NT);
% G. gradient operator: P1 -> (ND0)'
G = Me*icdmat(double(edge),[-1 1]);  % gradient matrix
% B = G';                              % -divergence matrix
% R. rotation operator
R = icdmat(double(elem2edge),[1 -1 1]);
C = R'*invMt*R;
A = [-Mv G'; G C];

%% Assemble right hand side
F = zeros(Ndof,1);
if ~isfield(pde,'f') || (isfield(pde,'f') && isreal(pde.f) && all(pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order is 3
end
if ~isempty(pde.f)
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ... 
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for k = 1:3
            i = locEdge(k,1); j = locEdge(k,2);
            % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
            phi_k = lambda(p,i)*Dlambda(:,:,j)-lambda(p,j)*Dlambda(:,:,i);
            rhs = dot(phi_k,fp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end
    bt = bt.*repmat(area,1,3);
    F = accumarray(N+elem2edge(:),bt(:),[Ndof 1]);
end
clear pxy fp bt rhs phi_k Dlambda

%% Boundary Conditions
[AD,F,u,sigma,freeNode,freeEdge,freeDof,isPureNeumann] = getbdHodgeLapE(F);
G  = G(freeEdge,freeNode);
C  = C(freeEdge,freeEdge);
Mv = Mv(freeNode,freeNode);
f = F(freeNode);
g = F(N+freeEdge);

%% Record assembling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the linear system
if isfield(option,'solver')
    method = upper(option.solver);
else
    method = 'DIAG';
end

% multigrid options 
option.mg.isFreeEdge = freeEdge; % needed in mg solver
option.mg.isFreeNode = freeNode; 
switch method
    case 'DIRECT'
        t = cputime;
        temp = zeros(Ndof,1);
        temp(freeDof) = A(freeDof,freeDof)\F(freeDof);
        sigma(freeNode) = temp(freeNode);
        u(freeEdge) = temp(freeEdge+N);
        residual = norm(F - AD*[sigma; u]);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'TRI' % GMRES with a triangular preconditioner 
        [sigma0,u0,info] =  tripreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
    case 'SDIAG' % GMRES with a signed diagonal preconditioner 
        [sigma0,u0,info] =  sdiagpreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
    case 'APPF' % GMRES with a signed diagonal preconditioner 
        [sigma0,u0,info] =  appfactpreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
    case 'DIAG' % MINRES with a diagonal preconditioner
        [sigma0,u0,info] =  diapreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
end
if isPureNeumann
   u = u - mean(u); % normalize u 
end

%% Output
eqn = struct('A',AD,'Mv',Mv,'Me',Me,'G',G,'C',C,'f',F,'edge',edge,...
             'freeNode',freeNode,'freeEdge',freeEdge);
info.assembleTime = assembleTime;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdHodgeLapE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,f,u,sigma,freeNode,freeEdge,freeDof,isPureNeumann] = getbdHodgeLapE(f)

    u = zeros(NE,1); sigma = zeros(N,1);    
    % Find boundary edges: Neumann
    isNeumann(elem2edge(bdFlag(:)==2)) = true;
    Neumann = edge(isNeumann,:); 
    
    % Neumann boundary condition: modify nodal dof values
    if ~isfield(pde,'gun') || (isnumeric(pde.gun) && all(pde.gun == 0))
        pde.gun = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.gun)
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
            gNp = pde.gun(ppxy);
            for igN = 1:2
                ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
            end
        end
        ge = ge.*repmat(el,1,2);
        f = f + accumarray(Neumann(:), ge(:),[Ndof,1]); 
    end
    % Neumann boundary condition: modify edge dof values
    edgeSign = ones(NE,1);
    idx = (bdFlag(:,1) ~= 0) & (elemSign == -1);% first edge is on boundary
    edgeSign(elem2edge(idx,1)) = -1;
    idx = (bdFlag(:,2) ~= 0) & (elemSign == 1); % second edge is on boundary
    edgeSign(elem2edge(idx,2)) = -1;
    idx = (bdFlag(:,3) ~= 0) & (elemSign == -1);% first edge is on boundary
    edgeSign(elem2edge(idx,3)) = -1;
    if ~isfield(pde,'grotu') || (isnumeric(pde.grotu) && all(pde.grotu==0))
        pde.grotu = [];
    end
    if ~isempty(pde.grotu) && ~isempty(Neumann)
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [lambda,weight] = quadpts1(option.gNquadorder);
        nQuad = size(lambda,1);
        Neumannidx = find(isNeumann) + N;
        for ip = 1:nQuad
        	pxy = lambda(ip,1)*node(Neumann(:,1),:)+...
                  lambda(ip,2)*node(Neumann(:,2),:);               
            f(Neumannidx) = f(Neumannidx) + weight(ip)*pde.grotu(pxy);
        end
        f(Neumannidx) = f(Neumannidx).*edgeSign(isNeumann);
        % no edge length since the basis of edge element contains it.
    end
    
    % Find Dirichlet boundary nodes and edges
    isDirichlet = false(NE,1);
    isDirichlet(elem2edge(bdFlag(:)==1)) = true;
    Dirichlet = edge(isDirichlet,:);
    isfixedNode = false(N,1); 
    isfixedNode(Dirichlet(:)) = true;
    fixedNode = find(isfixedNode);
    fixedDof = [fixedNode; N + find(isDirichlet)];
    freeNode = find(~isfixedNode);
    freeEdge = find(~isDirichlet);
    freeDof = [freeNode; N + freeEdge];
    isPureNeumann = false;
    if isempty(fixedNode) % pure Neumann boundary condition
        isPureNeumann = true;
        fixedDof = Ndof;
        freeDof = (1:Ndof-1)';    % eliminate the kernel
    end
    
    % Modify right hand side to include Dirichlet boundary condition
    if ~isfield(pde,'gu') || (isnumeric(pde.gu) && all(pde.gu == 0))   % zero gu
        pde.gu = [];
    end
    if ~isfield(pde,'gsigma') || (isnumeric(pde.gsigma) && all(pde.gsigma == 0))   % zero gsigma
        pde.gsigma = [];
    end
    if ~isPureNeumann && ~isempty(fixedNode) && (~isempty(pde.gu) || ~isempty(pde.gsigma))
        sigma = zeros(N,1);
        sigma(fixedNode) = pde.gsigma(node(fixedNode,:));
        u = zeros(NE,1);
        u(isDirichlet) = edgeinterpolate(pde.gu,node,Dirichlet);
        f = f - A*[sigma; u];
        f(fixedDof) = [sigma(fixedNode); u(isDirichlet)];    
    end
    
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0.
    if ~isempty(fixedDof)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedDof) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end 
    
    % Pure Neumann boundary condition
    if isPureNeumann
        f(N+1:Ndof) = f(N+1:Ndof) - mean(f(N+1:Ndof));   % compatilbe condition: sum(b) = 0
    end
    end
end
