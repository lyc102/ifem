function [u,edge,eqn,info] = eddycurrent(node,elem,bdFlag,pde,option)
%% EDDYCURRENT Eddy current equation: lowest order edge element.
%
% u = EDDYCURRENT(node,elem,bdFlag,pde) produces the lowest order edge
%   element approximation of the electric field of the following eddy
%   current equation.
%
% curl(mu^(-1)curl u) + beta u = J    in \Omega,  
%                        n × u = n × g_D  on \Gamma_D,
%          n × (mu^(-1)curl u) = n × g_N  on \Gamma_N.
% 
% which is from the time discretization of 
%        d_t(sigma E) + curl(mu^(-1)curl u) = - d_t j.
%
% The weak formulation is
% (mu^{-1}curl u, curl v) + (beta u,v) = (J,v) - <n × g_N,v>_{\Gamma_N}.
%
% It is the 2-D and positive definite version of Maxwell equation. For 3-D
% eddy current equation, use Maxwell.
%
% Reference: Beck, R. and Hiptmair, R. and Hoppe, R.H.W. and Wohlmuth, B.
% Mathematical Modelling and Numerical Analysis. 34(1):159--182, 2000.
%
% The data of the equation is enclosed in the pde structure:
%   - pde.mu      : permeability, i.e., magnetic constant/tensor
%   - pde.beta    : 
%   - pde.J       : current density
%   - pde.g_D     : Dirichlet boundary condition
%   - pde.g_N     : Neumann boundary condition
%
% The mesh is given by (node,elem). The boundary condition is specified by
% bdFlag; see <a href="matlab:ifem bddoc">bddoc</a>.
%
% The function eddycurrent assembes the matrix equation (A+M)*u = b and solves
% it by the direct solver (small size dof <= 2e3) or the HX preconditioned
% Krylov iterative methods (large size dof > 2e3).
% 
% u = EDDYCURRENT(node,elem,bdFlag,pde,option) specifies the solver options.
%   - option.solver == 'direct': the built in direct solver \ (mldivide)
%   - option.solver == 'mg':     multigrid-type solvers mg is used.
%   - option.solver == 'notsolve': the solution u = u_D. 
% The default setting is to use the direct solver for small size problems
% and multigrid solvers for large size problems. For more options on the
% multigrid solver mg, type help mg.
%
%
% Example
%
% See also Maxwell1, Maxwell2, cubeMaxwell, mgMaxwell
%
% Reference page in Help browser
%       <a href="matlab:ifem Maxwelldoc">Maxwelldoc</a> 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end

%% Construct Data Structure
[elem2dof,edge,dofSign] = dofedge(elem);
locEdge = [2,3;3,1;1,2];
N = size(node,1);   NT = size(elem,1);  Ndof = size(edge,1);

%% Compute coefficients
if ~isfield(pde,'mu'), pde.mu = 1; end
if ~isempty(pde.mu) && isnumeric(pde.mu)
    mu = pde.mu;                % mu is an array
else                            % mu is a function
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;
    mu = pde.mu(center);              
end
if ~isfield(pde,'beta'), pde.beta = 1; end
if ~isempty(pde.beta) && isnumeric(pde.beta)
    beta = pde.beta;      % beta is an array
else                            % beta is a function
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;
    beta = pde.beta(center);              
end

tstart = tic;
%% Element-wise basis
% edge indices of 3 local bases: 
% [2,3;3,1;1,2];
% phi = lambda_iDlambda_j - lambda_jDlambda_i;
% curl phi = 2*Dlambda_i × Dlambda_j;
[Dlambda,area] = gradbasis(node,elem);
curlPhi(:,3)=-2*(Dlambda(:,1,1).*Dlambda(:,2,2)-Dlambda(:,1,2).*Dlambda(:,2,1));
curlPhi(:,1)=-2*(Dlambda(:,1,2).*Dlambda(:,2,3)-Dlambda(:,1,3).*Dlambda(:,2,2));
curlPhi(:,2)=-2*(Dlambda(:,1,3).*Dlambda(:,2,1)-Dlambda(:,1,1).*Dlambda(:,2,3));
DiDj = zeros(NT,3,3);
for i = 1:3
    for j = i:3        
        DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2);
        DiDj(:,j,i) = DiDj(:,i,j);
    end
end

%% Assemble matrices
ii = zeros(6*NT,1); jj = zeros(6*NT,1); 
sA = zeros(6*NT,1); sM = zeros(6*NT,1);
index = 0;
for i = 1:3
    for j = i:3
        % local to global index map and its sign
        signij = dofSign(:,i).*dofSign(:,j);
        idx = (signij == -1); % inconsistent local direction
        % curl-curl matrix
        Aij = dot(curlPhi(:,i),curlPhi(:,j),2).*area./mu;
        Aij(idx) = -Aij(idx);
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));
        sA(index+1:index+NT) = Aij;
        % mass matrix
        % locEdge = [2,3;3,1;1,2];
        i1 = locEdge(i,1); i2 = locEdge(i,2);
        j1 = locEdge(j,1); j2 = locEdge(j,2);
        Mij = 1/12*area.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           - (1+(i1==j2))*DiDj(:,i2,j1) ...
                           - (1+(i2==j1))*DiDj(:,i1,j2) ...
                           + (1+(i2==j2))*DiDj(:,i1,j1));
        Mij = Mij.*beta;
        Mij(idx) = -Mij(idx);
        sM(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
clear curlPhi % clear large size data
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
M = M + MU + MU';
% bigA = A + M;

%% Assemble right hand side
f = zeros(Ndof,1);
if ~isfield(pde,'J') || (isfield(pde,'J') && isreal(pde.J) && all(pde.J==0))
    pde.J = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order is 3
end
if isfield(pde,'J') && ~isempty(pde.J)
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:);
        Jp = pde.J(pxyz);
        for k = 1:3
            i = locEdge(k,1); j = locEdge(k,2);
            % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
            phi_k = lambda(p,i)*Dlambda(:,:,j)-lambda(p,j)*Dlambda(:,:,i);
            rhs = double(dofSign(:,k)).*dot(phi_k,Jp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end
    bt = bt.*repmat(area,1,3);
    f = accumarray(elem2dof(:),bt(:),[Ndof 1]);
end
clear pxyz Jp bt rhs phi_k

%% Set up solver
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 1e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % Multigrid-type  solver for large size systems
        option.solver = 'cg';
    end
end
solver = option.solver;

%% Assembeling corresponding matrices for HX preconditioner
if ~strcmp(solver,'direct') && ~strcmp(solver,'nosolve')
    AP = sparse(N,N);  % AP = - div(mu^{-1}grad) + |Re(beta)| I
    BP = sparse(N,N);  % BP = - div(|Re(beta)|grad)
    for i = 1:3
        for j = i:3
            temp = DiDj(:,i,j).*area;
            Aij = 1./mu.*temp;
            Bij = abs(real(beta)).*temp;
            Mij = 1/12*abs(real(beta)).*area;
            if (j==i)
                AP = AP + sparse(elem(:,i),elem(:,j),Aij+2*Mij,N,N);
                BP = BP + sparse(elem(:,i),elem(:,j),Bij,N,N);            
            else
                AP = AP + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                                 [Aij+Mij; Aij+Mij],N,N);        
                BP = BP + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                                 [Bij; Bij],N,N);        
            end        
        end
    end
end
clear Aij Bij Mij

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
    bdFlag = setboundary(node,elem,'Dirichlet');
end
if ~isempty(bdFlag)
    % Find boundary edges and nodes
    isBdEdge = false(Ndof,1);
    isBdEdge(elem2dof(bdFlag(:,1) == 1,1)) = true;
    isBdEdge(elem2dof(bdFlag(:,2) == 1,2)) = true;
    isBdEdge(elem2dof(bdFlag(:,3) == 1,3)) = true;
    bdEdge = edge(isBdEdge,:);
    isBdNode(bdEdge) = true;
end
% modify the matrix to include the Dirichlet boundary condition
if any(isBdEdge)  % contains Dirichlet boundary condition
    bdidx = zeros(Ndof,1); 
    bdidx(isBdEdge) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    bigAD = T*(A+M)*T + Tbd;
    if ~strcmp(solver,'direct') && ~strcmp(solver,'nosolve')
        % modify the corresponding Poisson matrix
        bdidx = zeros(N,1); 
        bdidx(isBdNode) = 1;
        Tbd = spdiags(bdidx,0,N,N);
        T = spdiags(1-bdidx,0,N,N);
        AP = T*AP*T + Tbd;
        BP = T*BP*T + Tbd;
    end
else      % pure Neumann boundary condition
    bigAD = A + M;
    if ~strcmp(solver,'direct') && ~strcmp(solver,'nosolve')
       BP = BP + 1e-8*speye(N);  % make B non-singular      
    end
end

%% Part 2: Find boundary edges and modify the load b
g = zeros(Ndof,1);
% Find Neumann boundary faces
if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
    bdFlag = setboundary(node,elem,'Neumann');
end
% non-zero Neumann boundary condition
if ~isempty(bdFlag) && ~isempty(pde.g_N)
    % TODO: to be added
    f = f - g;
end
% nonzero Dirichlet boundary condition
u = zeros(Ndof,1);
if ~isempty(bdEdge) && ~isempty(pde.g_D) && ...
   ~(isnumeric(pde.g_D) && all(pde.g_D == 0))
    % else no bddof or g_D = 0 (no modification needed)
    if (isnumeric(pde.g_D) && length(pde.g_D) == Ndof)
        u(isBdEdge) = pde.g_D(isBdEdge);
    else
        u(isBdEdge) = edgeinterpolate(pde.g_D,node,bdEdge);
    end
    f = f - (A+M)*u;
    f(isBdEdge) = u(isBdEdge);
end
%% Remark
% The order of assign Neumann and Dirichlet boundary condition is
% important to get the right setting of the intersection of Dirichlet and
% Neumann faces.
    
%% Record assembling time
info.assembleTime = toc(tstart);
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 1
    fprintf('Time to assemble matrix equation %4.2g s\n',info.assembleTime);
end

%% Solve the system of linear equations
if strcmp(solver,'direct')
    % exact solver
    tstart = tic;
    freeDof = find(~isBdEdge);
    u(freeDof) = bigAD(freeDof,freeDof)\f(freeDof);
    time = toc(tstart); itStep = 0; flag = 2; err = norm(f - bigAD*u);    
elseif strcmp(solver,'nosolve')
    eqn = struct('A',A,'M',M,'f',f,'g',g,'bigA',bigAD,'isBdEdge',isBdEdge); 
    info = [];
    return;
else
%     u0 = edgeinterpolate(pde.g_D,node,edge);
    u0 = u;
    option.x0 = u0;
    [u,flag,itStep,err,time] = mgMaxwell(bigAD,f,AP,BP,node,elem,edge,[],isBdEdge,option);
end

%% Output
eqn = struct('A',A,'M',M,'f',f,'g',g,'bigA',bigAD,'isBdEdge',isBdEdge);
info = struct('solverTime',time,'itStep',itStep,'error',err,'flag',flag);