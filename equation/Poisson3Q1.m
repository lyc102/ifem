function [u,Du,eqn,info] = Poisson3Q1(node,elem,pde,bdFlag, option)
%% POISSON3Q1 Poisson equation: Q1 trilinear element.
%
%   u = Poisson3Q1(node,elem,pde,bdFlag) produces the trilinear finite element
%   approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
% Author: Huayi Wei < huayiwei1984@gmail.com>. 
%
% See also: PoissonQ1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


if ~exist('option','var'), option = []; end
[N,Dim] = size(node); 
[NT,NV] = size(elem);
Ndof = N;

tic;

hx = node(elem(:,2),1) - node(elem(:,1),1);
hy = node(elem(:,4),2) - node(elem(:,1),2);
hz = node(elem(:,5),3) - node(elem(:,1),3);
volume = abs(hx.*hy.*hz);

%% Assemble stiffness matrix
% generate sparse pattern
ii = zeros(36*NT,1); jj = zeros(36*NT,1); 
index = 0;
for i = 1:NV
    for j = i:NV
        ii(index+1:index+NT) = double(elem(:,i)); 
        jj(index+1:index+NT) = double(elem(:,j));  
        index = index + NT;
    end
end

% quadrature points
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder')
    option.dquadorder = 2;        % default order is exact for quadratic function
end
[GaussPt, weight] = quadptscube(option.dquadorder); % quadrature points on unit cube.
lambdax = GaussPt(:,1);
lambday = GaussPt(:,2);
lambdaz = GaussPt(:,3);

nQuad = size(GaussPt,1);
% compute non-zeros
sA = zeros(36*NT,nQuad);
for p = 1:nQuad
    if isfield(option,'isoparametric') && option.isoparametric
        [phi, Dphi, J] = hexbasis(node,elem,GaussPt(p,:));
    else
        Dphi = zeros(NT,3,8);
        phi2d = kron([lambday(p) 1-lambday(p)],[lambdax(p); 1-lambdax(p)]);
        phi = transpose(kron([lambdaz(p) 1-lambdaz(p)],phi2d([1 2 4 3])));
        Dx2d = kron([lambday(p) 1-lambday(p)],1./hx*[-1 1]);
        Dphi(:,1,:) = kron([lambdaz(p) 1-lambdaz(p)],Dx2d(:,[1 2 4 3]));
        Dy2d = kron(1./hy*[-1 1],[lambdax(p) 1-lambdax(p)]);
        Dphi(:,2,:) = kron([lambdaz(p) 1-lambdaz(p)],Dy2d(:,[1 2 4 3]));
        Dphi(:,3,:) = kron(1./hz*[-1 1],phi2d([1 2 4 3]));   
        J = volume;
    end
    index = 0;
    for i = 1:NV
        for j = i:NV
            Aij = 0;
            if isempty(pde.d) || isnumeric(pde.d)
                Aij = Aij + weight(p)*dot(Dphi(:,:,i),Dphi(:,:,j),2).*J;
            else
                pxyz = zeros(NT, Dim);
                for k = 1:Dim
                    xk = node(:,k);
                    pxyz(:,k) = xk(elem)*phi(:);
                end
                Aij = Aij + weight(p)*dot(Dphi(:,:,i),Dphi(:,:,j),2).*pde.d(pxyz).*J;
            end
            if ~isempty(pde.d) && isnumeric(pde.d) % d is piecewise constant
                Aij = pde.d.*Aij;
            end
            sA(index+1:index+NT,p) = Aij;
            index = index + NT;
        end
    end
end
sA = sum(sA,2);
% assemble the matrix
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj sA

%% Assemble the right hand side
b = zeros(N,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end

if ~isempty(pde.f) 
    [GaussPt, weight] = quadptscube(option.fquadorder);
    lambdax = GaussPt(:,1);
    lambday = GaussPt(:,2);
    lambdaz = GaussPt(:,3);    
    nQuad = size(GaussPt,1);
    bt = zeros(NT,NV);
    for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
        if isfield(option,'isoparametric') && option.isoparametric
            [phi, tempvar, J] = hexbasis(node,elem,GaussPt(p,:));
        else
            phi2d = kron([lambday(p) 1-lambday(p)],[lambdax(p); 1-lambdax(p)]);
            phi = transpose(kron([lambdaz(p) 1-lambdaz(p)],phi2d([1 2 4 3])));
            J = volume;
        end
%         [phi, tempvar, J] = hexbasis(node,elem, pts(p,:));
        pxyz = zeros(NT, Dim);
        for k = 1:Dim
            xk = node(:,k);
            pxyz(:,k) = xk(elem)*phi(:);
        end
		fp = pde.f(pxyz);
        bt = bt + (weight(p)*fp.*J)*phi';
    end
    b = accumarray(elem(:),bt(:),[Ndof 1]);
end
clear pxyz bt

%% Set up boundary conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,b,u,freeNode,isPureNeumann] = getbd3(b);

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
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
%         option.x0 = u;
        option.solver = 'CG';
%         [u,info] = mg(AD,b,elem,option);
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option);                 
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    patchArea = accumarray(elem(:),[J;J;J;J;J;J;J;J]/8, [N 1]); % TODO: Here is J?
    uc = sum(u.*patchArea)/sum(J);
    u = u - uc;   % int u = 0
end

%% Output information
eqn = struct('A',AD,'b',b,'freeNode',freeNode);
info.assembleTime = assembleTime;

%% Compute Du
Du = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbd3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeNode,isPureNeumann] = getbd3(b)
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
        allFace = [elem(:,[1 4 3 2]); elem(:,[1 2 6 5]); elem(:,[5 6 7 8]);...
                   elem(:,[8 7 3 4]); elem(:,[4 1 5 8]); elem(:,[2 3 7 6])];
        Robin = allFace(isRobin,:);
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        %% To do: this is for tet mesh. Modify for hex mesh later.
        v12 = node(Robin(:,2),:)-node(Robin(:,1),:);
        v13 = node(Robin(:,3),:)-node(Robin(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        if ~isfield(option,'gRquadorder')
            option.gRquadorder = 3;   % default order exact for linear gR
        end
        [lambdagR,weightgR] = quadpts(option.gRquadorder);
        nQuadgR = size(lambdagR,1);
        bdphi = lambdagR;
        NR = size(Robin,1);
        ss = zeros(NR,3,3);
        % int g_R phi_i phi_j
        for pp = 1:nQuadgR
            ppxyz = lambdagR(pp,1)*node(Robin(:,1),:) ...
                  + lambdagR(pp,2)*node(Robin(:,2),:) ...
                  + lambdagR(pp,3)*node(Robin(:,3),:);
            gRp = pde.g_R(ppxyz);
            for iR = 1:3
                for jR = iR:3
                    ss(:,iR,jR) = ss(:,iR,jR) + ...
                    weightgR(pp)*gRp*bdphi(pp,iR).*bdphi(pp,jR);
                end
            end
        end
        ss(:) = ss(:).*repmat(area,9,1);       
        index = 0;
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

    % Find Dirichlet boundary nodes: fixedNode
    fixedNode = []; freeNode = [];
    if ~isempty(bdFlag) % bdFlag specifies different bd conditions
        [fixedNode,bdFace,isBdNode] = findboundary3(elem,bdFlag);
        freeNode = find(~isBdNode);
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given
        [fixedNode,bdFace,isBdNode] = findboundary3(elem);
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
    bdidx = zeros(Ndof,1); 
    bdidx(fixedNode) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    AD = T*A*T + Tbd;

    %% Part 2: Find boundary faces and modify the load b
    % Find boundary faces: Neumann
    Neumann = [];
    if ~isempty(bdFlag)  % bdFlag specifies different bd conditions
        Neumann = bdFace;       
    end
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        [tempvar,Neumann] = findboundary3(elem); %#ok<ASGLU>
    end

    % Neumann boundary condition
    if ~isempty(Neumann) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
        v14 = node(Neumann(:,4),:)-node(Neumann(:,1),:);
        area = sqrt(abs(sum(mycross(v12,v14,2).^2,2)));
        % three middle points rule
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [GaussPt,weightgN] = quadptsquad(option.gNquadorder);
        xxi = GaussPt(:,1);
        eta = GaussPt(:,2);
        nQuadgN = size(GaussPt,1);
        % barycentrical coordinate for four basis
        lambdagN = zeros(nQuadgN,4);
        lambdagN(:,1) = (1 - xxi).*(1 - eta);
        lambdagN(:,2) = xxi.*(1 - eta);
        lambdagN(:,3) = xxi.*eta;
        lambdagN(:,4) = (1 - xxi).*eta;
        ge = zeros(size(Neumann,1),4);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxyz = node(Neumann(:,1),:)*lambdagN(pp,1) + node(Neumann(:,2),:)*lambdagN(pp,2)...
                  + node(Neumann(:,3),:)*lambdagN(pp,3) + node(Neumann(:,4),:)*lambdagN(pp,4);
            gNp = pde.g_N(ppxyz);
            for iN = 1:4
                ge(:,iN) = ge(:,iN) + weightgN(pp)*lambdagN(pp,iN)*gNp;
            end
        end
        ge = ge.*repmat(area,1,4);
        b = b + accumarray(Neumann(:),ge(:),[Ndof,1]); 
    end
    % The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann edges and no modification of
    % A,u,b is needed.

    % Dirichlet boundary condition
    if ~isPureNeumann && ~isempty(fixedNode) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
        if isnumeric(pde.g_D)
            u(fixedNode) = pde.g_D(fixedNode);
        else
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
    end % end of getbd3
end % end of Poisson3Q1
