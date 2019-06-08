function [u,Du,eqn,info] = PoissonQ1old(node,elem,pde,bdFlag, option)
%% POISSONQ1 Poisson equation: Q1 bilinear element.
%
%   u = PoissonQ1(node,elem,pde,bdFlag) produces the bilinear finite element
%   approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
%   The quad mesh is given by node and elem and the boundary edge is given by
%   bdFlag. See meshdoc, bddoc for details. The data is given by the
%   structure pde which contains function handles f, g_D, g_N, g_R, or d.
%   For general elliptic equations with convection and reaction
%   coefficients, see ellipticpde.
%   
%   The function Poisson assembes the matrix equation AD*u = b and solves
%   it by the direct solver.
%
%   The diffusion coefficient d is a scalar function or a column array with
%   the same length as the elem array. 
% 
%   u = PoissonQ1(node,elem,pde,bdFlag,option) specifies the options.
%    - option.dquadorder: quadrature order for diffusion coefficients
%    - option.fquadorder: quadrature order for computing right hand side f
%
%   When only one type of boundary condition is imposed, the input argument
%   bdFlag can be skipped. The boundary condition is implicitly given in
%   the pde structure by specifying g_D or g_N only. See examples below.
%
%   [u,Du] = PoissonQ1(node,elem,pde,bdFlag) returns also the
%   non-modified stiffness matrix A, which is semi-definite. The kernel of
%   A consists of constant vectors. The matrix A can be used to evulate the
%   bilinear form A(u,v) = u'*A*v, especially the enery norm of a finite
%   element function u by sqrt(u'*A*u).
%
%   [u,Du,eqn] = PoissonQ1(node,elem,pde,bdFlag) returns also the equation
%   structure eqn, which includes: 
%     - eqn.AD:  modified stiffness matrix AD;
%     - eqn.b:   the right hand side. 
%   The solution u = AD\b. The output eqn can be used to test other solvers.
%
%
%   Example
%     clear all
%     node=  [0, 0; 1, 0; 1, 1; 0,1];
%     elem = [1,2,3,4];
%     for k = 1:4
%         [node,elem] = uniformrefinequad(node,elem);
%     end
%     % Homogenous Dirichlet boundary condition
%     pde.f = inline('ones(size(p,1),1)','p');
%     pde.g_D = 0;
%     u = PoissonQ1(node,elem,pde);
%     figure(1); 
%     showsolution(node,elem,u);
%     % Non-homogenous Dirichlet boundary condition
%     pde.f = inline('-4*ones(size(p,1),1)','p');
%     pde.g_D = inline('p(:,1).^2 + p(:,2).^2','p');
%     u = PoissonQ1(node,elem,pde);
%     figure(2); 
%     showsolution(node,elem,u);
%     % Homogenous Neumann boundary condition
%     clear pde
%     pde.f = inline('pi^2*cos(pi*p(:,1)).*cos(pi*p(:,2))','p');
%     u = PoissonQ1(node,elem,pde);
%     figure(3);
%     showsolution(node,elem,u);
%
% See also: Poisson3Q1
% 
% Author: Huayi Wei < huayiwei1984@gmail.com>. 
%
% Modified by Long Chen. Change the assembling process. Add mg solver and
% update boundary condition.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
[N,Dim] = size(node); 
[NT,NV] = size(elem);
Ndof = N;

tic;
%% Assemble stiffness matrix
% generate sparse pattern
ii = zeros(10*NT,1); jj = zeros(10*NT,1); 
index = 0;
for i = 1:4
    for j = i:4
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
[pts, w] = quadptsquad(option.dquadorder);
nQuad = size(pts,1);
% compute non-zeros
sA = zeros(10*NT,nQuad);
for p = 1:nQuad
    % Dphi at quadrature points
   [phi, Dphip, J] = quadbasis(node,elem,pts(p,:)); 
    index = 0;
    for i = 1:4
        for j = i:4
            Aij = 0;
            if isempty(pde.d) || isnumeric(pde.d)
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            else
                pxy = zeros(NT, Dim);
                for ip = 1:Dim
                    xi = node(:,ip);
                    pxy(:,ip) = xi(elem)*phi;
                end
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.d(pxy);
            end
            if ~isempty(pde.d) && isnumeric(pde.d) % d is piecewise constant
                Aij = pde.d.*Aij;
            end
            Aij = Aij.*J;
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
clear Aij ii jj Dphip

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end

if ~isempty(pde.f) 
    [pts, weight] = quadptsquad(option.fquadorder);
    nQuad = size(pts,1);
    bt = zeros(NT,NV);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
        [phi, tempvar, J] = quadbasis(node,elem, pts(p,:));
        pxy = zeros(NT, Dim);
        for i = 1:Dim
            xi = node(:,i);
            pxy(:,i) = xi(elem)*phi; % ? questionable 
        end
		fp = pde.f(pxy);
        bt = bt + (weight(p)*fp.*J)*phi';
    end
    b = accumarray(elem(:),bt(:),[N 1]);
end
clear pxy bt

%% Set up boundary conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,b,u,freeNode,isPureNeumann] = getbd(b);

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
        [u,info] = mg(AD,b,elem,option);
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option);                 
end
% post-process for pure Neumann problem
if isPureNeumann
    patchArea = accumarray(elem(:),[J;J;J;J]/4, [N 1]); 
    uc = sum(u.*patchArea)/sum(J);
    u = u - uc;   % int u = 0
end

%% Output information
eqn = struct('A',AD,'b',b,'freeNode',freeNode);
info.assembleTime = assembleTime;

%% Compute Du
Du = [];

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
        allEdge = [elem(:,[1,2]); elem(:,[2,3]); elem(:,[3,4]); elem(:,[4 1])];    % quad
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
%         ss = [1/4*temp, 1/4*temp, 1/4*temp, 1/4*temp];
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
    end
    end % end of getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % end of PoissonQ1