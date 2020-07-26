function [soln,eqn,info] = Maxwell1saddle(node,elem,bdFlag,pde,option,varargin)
%% MAXWELL1SADDLE Maxwell equation in mixed form: linear edge element.
%
% u = Maxwell1saddle(node,elem,bdFlag,pde) produces the lowest order edge
%   element approximation of the electric field of the magnetostatics
%
%                     curl(mu^(-1)curl u)  = J    in \Omega,  
%                                  div  u  = g    in \Omega,  
%                                   n x u = n x g_D  on \Gamma_D,
%                     n x (mu^(-1)curl u) = n x g_N  on \Gamma_N.
% 
% based on the weak formulation
%
% (mu^{-1}curl u, curl v) +(grad p ,v)= (J,v) - <n x g_N,v>_{\Gamma_N}.
% (u,grad q)                          = -(g,q)                    
%
% 
% The data of the equation is enclosed in the pde structure:
%   - pde.mu      : permeability, i.e., magnetic constant/tensor
%   - pde.J       : current density
%   - pde.g       : divergence
%   - pde.g_D     : Dirichlet boundary condition
%   - pde.g_N     : Neumann boundary condition
%
% 
% u = Maxwell1saddle(node,elem,bdFlag,pde,option) specifies the solver options.
%   - option.solver == 'direct': the built in direct solver \ (mldivide)
%   - option.solver == 'mg':     multigrid-type solvers mg is used.
%   - option.solver == 'none':   the solution u = u_D. 
%
% The default setting is to use the direct solver for small size problems
% and multigrid solvers for large size problems. For more options on the
% multigrid solver mg, type help mg.
%
% [soln,eqn] = Maxwell1saddle(node,elem,bdFlag,pde) returns soln structure
% which contains soln.u and soln.p and the equation structure eqn: 
% - eqn.A: matrix for differential operator;
% - eqn.M: mass matrix;
% - eqn.f: right hand side 
% - eqn.g: vector enclosed the Neumann boundary condition
% - eqn.edge: edge array
% - eqn.isFreeEdge: not Dirichlet boundary edges
%
% [soln,eqn] = Maxwell1saddle(node,elem,pde,bdEdge) returns also the
% information on the assembeling and solver, which includes:
% - info.assembleTime: time to assemble the matrix equation
% - info.solverTime:   time to solve the matrix equation
% - info.itStep:       number of iteration steps for the mg solver
% - info.error:        l2 norm of the residual b - A*u
% - info.flag:         flag for the mg solver.
%   flag = 0: converge within max iteration 
%   flag = 1: iterated maxIt times but did not converge
%   flag = 2: direct solver
%   flag = 3: no solve
%
% Example
%
%   cubeMaxwellSaddle1
%
% See also Maxwellsaddle, mgMaxwellsaddle
%
% Reference page in Help browser
%       <a href="matlab:ifem MaxwellsaddleND0femrate">MaxwellsaddleND0femrate</a> 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Set up optional input arguments
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
t = cputime;

%% Sort elem to ascend ordering
elemMG  = elem;     % save elem and bdFlag for multigrid
bdFlagMG = bdFlag;
[elem,bdFlag] = sortelem3(elem,bdFlag);

%% Construct data structure
[elem2edge,edge] = dof3edge(elem);
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
N = size(node,1);   NT = size(elem,1);  NE = size(edge,1);   
Nudof = 2*NE;    Npdof = N + NE;
elem2dofP2 = [elem N+elem2edge];

%% Compute coefficients
if ~isfield(pde,'mu'), pde.mu = 1; end
if ~isempty(pde.mu) && isnumeric(pde.mu)
    mu = pde.mu;                % mu is an array
else                            % mu is a function
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;
    mu = pde.mu(center);              
end
if ~isfield(pde,'epsilon')
    pde.epsilon = []; 
    epsilon = 0;
end
if isfield(pde,'omega')
    omega = pde.omega;
else
    omega = 1;
end
if ~isempty(pde.epsilon) 
    if isnumeric(pde.epsilon)
        epsilon = pde.epsilon;      % epsilon is an array
    else                            % epsilon is a function
        center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
                  node(elem(:,3),:) + node(elem(:,4),:))/4;
        epsilon = pde.epsilon(center);              
    end
    epsilon = omega^2*epsilon; 
end

%% Element-wise basis
% edge indices of 6 local bases: 
% [1 2], [1 3], [1 4], [2 3], [2 4], [3 4]
% phi = lambda_i Dlambda_j - lambda_j Dlambda_i;
% curl phi = 2*Dlambda_i \times Dlambda_j;
% psi = lambda_iDlambda_j + lambda_jDlambda_i;
% curl psi = 0;

[Dlambda,volume] = gradbasis3(node,elem);
curlPhi(:,:,6) = 2*mycross(Dlambda(:,:,3),Dlambda(:,:,4),2);
curlPhi(:,:,1) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,2),2);
curlPhi(:,:,2) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,3),2);
curlPhi(:,:,3) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,4),2);
curlPhi(:,:,4) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,3),2);
curlPhi(:,:,5) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,4),2);
DiDj = zeros(NT,4,4);
for i = 1:4
    for j = i:4        
        DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2);
        DiDj(:,j,i) = DiDj(:,i,j);
    end
end

%% Assemble matrices
ii = zeros(21*NT,1); jj = zeros(21*NT,1); 
sA = zeros(21*NT,1); sMphi = zeros(21*NT,1); sMpsi = zeros(21*NT,1);
index = 0;
for i = 1:6
    for j = i:6
        % curl-curl matrix % no contribution from psi
        Aij = dot(curlPhi(:,:,i),curlPhi(:,:,j),2).*volume./mu;
        ii(index+1:index+NT) = double(elem2edge(:,i)); 
        jj(index+1:index+NT) = double(elem2edge(:,j));
        sA(index+1:index+NT) = Aij;
        % mass matrix
        % locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        i1 = locEdge(i,1); i2 = locEdge(i,2);
        j1 = locEdge(j,1); j2 = locEdge(j,2);
        % (phi_i,phi_j)
        Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           - (1+(i1==j2))*DiDj(:,i2,j1) ...
                           - (1+(i2==j1))*DiDj(:,i1,j2) ...
                           + (1+(i2==j2))*DiDj(:,i1,j1));
        sMphi(index+1:index+NT) = Mij;
        % (psi_i,psi_j)
        Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           + (1+(i1==j2))*DiDj(:,i2,j1) ...
                           + (1+(i2==j1))*DiDj(:,i1,j2) ...
                           + (1+(i2==j2))*DiDj(:,i1,j1));
        sMpsi(index+1:index+NT) = Mij;
        % increment of index               
        index = index + NT;
    end
end
clear curlPhi % clear large size data
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Nudof,Nudof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Nudof,Nudof);
A = A + AU + AU';
M = sparse([ii(diagIdx); ii(diagIdx)+NE], [jj(diagIdx); jj(diagIdx)+NE],...
           [sMphi(diagIdx); sMpsi(diagIdx)],Nudof,Nudof);
MU = sparse([ii(upperIdx); ii(upperIdx)+NE], [jj(upperIdx); jj(upperIdx)+NE],...
           [sMphi(upperIdx); sMpsi(upperIdx)],Nudof,Nudof);
M = M + MU + MU';
clear AU MU
% (psi_i,phi_j)
ii = zeros(36*NT,1); jj = zeros(36*NT,1); ss = zeros(36*NT,1);
index = 0;
for i = 1:6
    for j = 1:6
        % local to global index map and its sign
        i1 = locEdge(i,1); i2 = locEdge(i,2);
        j1 = locEdge(j,1); j2 = locEdge(j,2);
        Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           - (1+(i1==j2))*DiDj(:,i2,j1) ...
                           + (1+(i2==j1))*DiDj(:,i1,j2) ...
                           - (1+(i2==j2))*DiDj(:,i1,j1));
%         Mij = Mij.*epsilon;
        ii(index+1:index+NT) = double(elem2edge(:,i))+NE; 
        jj(index+1:index+NT) = double(elem2edge(:,j));
        ss(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
ML = sparse(ii,jj,ss,Nudof,Nudof);
M  = M  + ML + ML';
clear ML
if abs(epsilon) > eps
    A = A - epsilon*M;
end

%% grad matrix       
grad = icdmat(double(edge),[-1,1]); % P1 to ND0
% e2v matrix: P1 to P2 lambda_i -> lambda_i(2lambda_i - 1)
ii = [(1:NE)';(1:NE)'];
jj = double(edge(:));
ss = -2*[ones(NE,1),ones(NE,1)];
e2v = sparse(ii,jj,ss,NE,N);
bigGrad = [grad, sparse(NE,NE); ...
            e2v, 4*speye(NE,NE)];
G = M*bigGrad;

%% Assemble right hand side
f = zeros(Nudof,1);
if ~isfield(pde,'J') || (isfield(pde,'J') && isreal(pde.J) && all(pde.J==0))
    pde.J = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order is 3
end
if isfield(pde,'J') && ~isempty(pde.J) && ~isnumeric(pde.J)
    [lambda,w] = quadpts3(option.fquadorder);
    nQuad = size(lambda,1);
    bt = zeros(NT,12);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:) ... 
             + lambda(p,4)*node(elem(:,4),:);
        Jp = pde.J(pxyz);
    %   locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        for k = 1:6
            i = locEdge(k,1); j = locEdge(k,2);
            % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
            phi_k = lambda(p,i)*Dlambda(:,:,j)-lambda(p,j)*Dlambda(:,:,i);
            rhs = dot(phi_k,Jp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
            % psi_k = lambda_iDlambda_j + lambda_jDlambda_i;
            psi_k = lambda(p,i)*Dlambda(:,:,j)+lambda(p,j)*Dlambda(:,:,i);
            rhs = dot(psi_k,Jp,2);
            bt(:,k+6) = bt(:,k+6) + w(p)*rhs;            
        end
    end
    bt = bt.*repmat(volume,1,12);
    f = accumarray([elem2edge(:); elem2edge(:)+NE],bt(:),[Nudof 1]);
% elseif isfield(pde,'J') && ~isempty(pde.J) && isnumeric(pde.J)
%     switch size(pde.J,1)
%         case NE % rhs already computed
%             f = pde.J;
%         case NT % piecwise constant
%             bt = zeros(NT,6);
%             for k = 1:6
%                 i = locEdge(k,1); j = locEdge(k,2);
%                 % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
%                 phi_k = (Dlambda(:,:,j)-Dlambda(:,:,i))/4;
%                 bt(:,k) = dot(phi_k,pde.J,2);
%             end
%             bt = bt.*repmat(volume,1,6);
%             f = accumarray(elem2dof(:),bt(:),[NE 1]);
%     end
end
clear pxyz Jp bt rhs phi_k psi_k

%% if u is not divergence free g0 = div u
g0 = zeros(Npdof,1);
if isfield(pde,'g') && ~isempty(pde.g) && ~isnumeric(pde.g)
    [lambda,w] = quadpts3(option.fquadorder+1);
    nQuad = size(lambda,1);
    gt = zeros(NT,10);
    phi(:,10)= 4*lambda(:,3).*lambda(:,4);
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phi(:,4) = lambda(:,4).*(2*lambda(:,4)-1);
    phi(:,5) = 4*lambda(:,1).*lambda(:,2);
    phi(:,6) = 4*lambda(:,1).*lambda(:,3);
    phi(:,7) = 4*lambda(:,1).*lambda(:,4);
    phi(:,8) = 4*lambda(:,2).*lambda(:,3);
    phi(:,9) = 4*lambda(:,2).*lambda(:,4);
    for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
		gp = pde.g(pxyz);
        for j = 1:10
            gt(:,j) = gt(:,j) + w(p)*phi(p,j)*gp;
        end
    end
    gt = gt.*repmat(volume,1,10);
    g0 = -accumarray(elem2dofP2(:),gt(:),[Npdof 1]); % - div u = g0
end

%% Record assembling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Boundary conditions
[AD,f,u,isFreeNode,isFreeEdge] = getbdMaxwell1saddle(f);
p  = zeros(Npdof,1);
g0 = g0 - G'*u;
% find free dofs for u and p
isFreeDofu = [isFreeEdge; isFreeEdge];
isFreeDofp = [isFreeNode; isFreeEdge];

%% Output
eqn = struct('A',A,'M',M,'G',G,'f',f,'g',g0,'AD',AD,'edge',edge,'isFreeEdge',isFreeEdge);
info.assembleTime = assembleTime;

%% reduce to free dofs
A  = A(isFreeDofu,isFreeDofu);
G  = G(isFreeDofu,isFreeDofp);
bigGrad = bigGrad(isFreeDofu,isFreeDofp);
M = M(isFreeDofu,isFreeDofu);
f = f(isFreeDofu);
g0 = g0(isFreeDofp);

%% Solve the linear system
if isfield(option,'solver')
    method = upper(option.solver);
else
    method = 'MG';
end
if nargin>=6
    HB = varargin{1};
else
    HB = [];
end        

% multigrid options 
option.mg.isFreeEdge = isFreeEdge; % needed in mg solver
option.mg.isFreeNode = isFreeNode; 
switch method
    case 'DIRECT'
        t = cputime;
        Npi = size(G,2); 
        Nui = size(G,1);
        bigA   = [A  G; ...
                  G' sparse(Npi,Npi)];
        temp   = bigA\[f; g0];
        u(isFreeDofu) = temp(1:Nui);
        p(isFreeDofp) = temp(Nui+1:Nui+Npi);
        residual = norm([f;g0] - bigA*temp);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'MG'
        [u0,p0,info] = mgMaxwellsaddle(A,G,f,g0,node,elemMG,bdFlagMG,M,grad,option.mg,HB);
        u(isFreeDofu)  = u0;
        p(isFreeDofu)  = p0;        
end

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'p',p);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdMaxwell1saddle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,f,u,isFreeNode,isFreeEdge] = getbdMaxwell1saddle(f)

    %% Boundary conditions
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end
    if (isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R))
        % no boundary data is given = homogenous Neumann boundary condition
        bdFlag = []; 
    end

    %% Part 1: Find Dirichlet dof and modify the matrix
    % Find Dirichlet boundary dof: fixedDof
    isBdEdge = [];
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N)
        % Dirichlet boundary condition only
        bdFlag = setboundary3(node,elem,'Dirichlet');
    end
    if ~isempty(bdFlag)
        % Find boundary edges and nodes
        isBdEdge = false(NE,1);
        isBdNode = false(N,1);
        isBdEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
        isBdEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
        isBdEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
        isBdEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
        bdDof = [find(isBdEdge); NE + find(isBdEdge)];
        bdEdge = edge(isBdEdge,:);
        isBdNode(bdEdge) = true;
    end
    % modify the matrix to include the Dirichlet boundary condition
    if any(isBdEdge)  % contains Dirichlet boundary condition
        bdidx = zeros(Nudof,1); 
        bdidx(bdDof) = 1;
        Tbd = spdiags(bdidx,0,Nudof,Nudof);
        Te = spdiags(1-bdidx,0,Nudof,Nudof);
        AD = Te*A*Te + Tbd;
    else      % pure Neumann boundary condition
        isPureNeumann = 1;
        AD = A;
    end

    %% Part 2: Find boundary edges and modify the load b
    g = zeros(Nudof,1);
    % Find Neumann boundary faces
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        bdFlag = setboundary3(node,elem,'Neumann');
    end
    % non-zero Neumann boundary condition
    if ~isempty(bdFlag) && ~isempty(pde.g_N)
        % face 1
        isBdElem = find(bdFlag(:,1) == 2); %#ok<*NASGU>
        face = [2 3 4]; face2locdof = [6 5 4];
        if ~isempty(isBdElem)
            bdb = bdfaceintegral1(isBdElem,face,face2locdof);
            g = bdb; 
        end
        % face 2
        isBdElem = find(bdFlag(:,2) == 2);
        face = [1 3 4]; face2locdof = [6 3 2];
        if ~isempty(isBdElem)
            bdb = bdfaceintegral1(isBdElem,face,face2locdof);
            g = g + bdb;
        end
        % face 3
        isBdElem = find(bdFlag(:,3) == 2);
        face = [1 2 4]; face2locdof = [5 3 1];
        if ~isempty(isBdElem)
            bdb = bdfaceintegral1(isBdElem,face,face2locdof);
            g = g + bdb;
        end
        % face 4
        isBdElem = find(bdFlag(:,4) == 2);
        face = [1 2 3]; face2locdof = [4 2 1];
        if ~isempty(isBdElem)
            bdb = bdfaceintegral1(isBdElem,face,face2locdof);
            g = g + bdb;
        end
        f = f - g;
    end
    % nonzero Dirichlet boundary condition
    u = zeros(Nudof,1);
    if ~isempty(bdEdge) && ~isempty(pde.g_D) && ...
       ~(isnumeric(pde.g_D) && all(pde.g_D == 0))
        % else no bddof or g_D = 0 (no modification needed)
        if (isnumeric(pde.g_D) && length(pde.g_D) == Nudof)
            u(bdDof) = pde.g_D(bdDof);
        else
            u(bdDof) = edgeinterpolate1(pde.g_D,node,bdEdge);
        end
        f = f - A*u;
        f(bdDof) = u(bdDof);
    end
    % find free node, free edge and free dof
    isFreeEdge = ~isBdEdge;
    isFreeNode = ~isBdNode;
    %% Remark
    % The order of assign Neumann and Dirichlet boundary condition is
    % important to get the right setting of the intersection of Dirichlet and
    % Neumann faces.
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions bdfaceintegral1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function bdb = bdfaceintegral1(isBdElem,face,face2locdof)
    %% Compute boundary surface integral of lowest order edge element.
    %  bdb(k) = \int_{face} (n×g_N, phi_k) dS
    %  bdb(k+3) = \int_{face} (n×g_N, psi_k) dS
    
    %% Compute scaled normal
    faceIdx = true(4,1);
    faceIdx(face) = false;
    normal = -3*repmat(volume(isBdElem),1,3).*Dlambda(isBdElem,:,faceIdx);
    
    %% Data structure
    tetLocEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % edge of a tetrahedral [1 2 3 4]
    face2locEdge = [2 3; 1 3; 1 2]; % edge of the face [1 2 3]

    %% Compute surface integral
    Nbd = length(isBdElem);
    bt = zeros(Nbd,6);
    idx = zeros(Nbd,6,'int32');
    [lambda,w] = quadpts(3); % quadrature order is 3
    nQuad = size(lambda,1);
    for pp = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(pp,1)*node(elem(isBdElem,face(1)),:) ...
             + lambda(pp,2)*node(elem(isBdElem,face(2)),:) ... 
             + lambda(pp,3)*node(elem(isBdElem,face(3)),:);
        gNp = pde.g_N(pxyz,normal);    
        for s = 1:3
            kk = face2locdof(s);
            phi_k = lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                  - lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,1));
            rhs = dot(phi_k,gNp,2);
            bt(:,s) = bt(:,s) + w(pp)*rhs; % area is included in normal; see line 28
            idx(:,s) = elem2edge(isBdElem,kk);
            % psi_k = lambda_iDlambda_j + lambda_jDlambda_i;
            psi_k = lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                  + lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,1));
            rhs = dot(psi_k,gNp,2);
            bt(:,s+3) = bt(:,s+3) + w(pp)*rhs; % area is included in normal; see line 27
            idx(:,s+3) = elem2edge(isBdElem,kk) + NE;
        end
    end
    %% Distribute to DOF
    bdb = accumarray(idx(:),bt(:),[Ndof 1]);        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end