function [u,edge,eqn,info] = Maxwell1H(node,elem,HB,pde,bdFlag,option)
%% MAXWELL1H Maxwell equation for H: linear edge element.
%
% u = Maxwell1H(node,elem,HB,pde,bdFlag) produces the lowest order edge
%   element approximation of the magnetic field H of the time harmonic
%   Maxwell equation.
%
% curl(1/epsilon curl H) - omega^2*mu H = curl(1/epsilon J)    in \Omega,  
%                                 n × H = n × g_D  on \Gamma_D,
%    n × (1/epsilon curl H) = n×(-i*omega*E + 1/epsilon J)  on \Gamma_N.
% 
% based on the weak formulation
%
%   (1/epsilon curl u, curl v) - (omega^2*mu u,v) 
% = (1/epsilon J,curl v) - <n × -i*omega*E,v>_{\Gamma_N}.
%
% Here we use the fact curl H = -i*omega*epsilon*E + J. The boundary trace
% for J is canceled and the Neumman data g_N = curl H /epsilon = -i*omega*E.
% 
% The data is enclosed in the structure pde
%   - pde.mu      : permeability, i.e., magnetic constant
%   - pde.epsilon : a complex dielectric constant
%   - pde.omega   : wave number
%   - pde.J       : current density
%   - pde.g_D     : Dirichlet boundary condition
%   - pde.g_N     : Neumann boundary condition
%
% The mesh is given by (node,elem) and HB is needed for fast solver. The
% boundary faces is specified by bdFlag; see <a href="matlab:ifem bddoc">Data Structure: Boundary conditions.</a>  
%
% The function Maxwell1 assembes the matrix equation (A+M)*u = b and solves
% it by the direct solver (small size < 2e3) or the HX preconditioned CG
% iterative method (large size >=2e3).
% 
% u = Maxwell1(node,elem,HB,pde,bdFlag,option) specifies the solver options.
%   - option.solver == 'direct': the built in direct solver \ (mldivide)
%   - option.solver == 'mg':     multigrid-type solvers mg is used.
%   - option.solver == 'notsolve': the solution u = u_D. 
% The default setting is to use the direct solver for small size problems
% and multigrid solvers for large size problems. For more options on the
% multigrid solver mg, type help mg.
%
% [u,edge] = Poisson(node,elem,pde,bdEdge) returns also the edge array
% which is essential for edge elements. 
%
% [u,edge,eqn] = Poisson(node,elem,pde,bdEdge) returns also the equation
% structure eqn, which includes: 
% - eqn.A: matrix for differential operator;
% - eqn.M: mass matrix;
% - eqn.f: right hand side 
% - eqn.g: vector enclosed the Neumann boundary condition
%
% Example
%   cubeMaxwell1
%
% See also Maxwell, Maxwell2, cubeMaxwell1, mgMaxwell
%
% Reference page in Help browser
%       <a href="matlab:ifem Maxwelldoc">doc Maxwell1</a> 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
tic;

%% Construct Data Structure
[elem2dof,dofSign,edge] = dof3edge(elem);
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
N = size(node,1);   NT = size(elem,1);  NE = size(edge,1);  Ndof = 2*NE;

%% Compute coefficients
if ~isfield(pde,'mu'), pde.mu = 1; end
if ~isempty(pde.mu) && isnumeric(pde.mu)
    mu = pde.mu;                % mu is an array
else                            % mu is a function
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;
    mu = pde.mu(center);              
end
if ~isfield(pde,'epsilon'), pde.epsilon = 1; end
if ~isempty(pde.epsilon) && isnumeric(pde.epsilon)
    epsilon = pde.epsilon;      % epsilon is an array
else                            % epsilon is a function
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;
    epsilon = pde.epsilon(center);              
end
if isfield(pde,'omega')
    omega = pde.omega;
    mu = omega^2*mu; 
else
    omega = 1;
end

%% Element-wise basis
% [1 2], [1 3], [1 4], [2 3], [2 4], [3 4]
% phi = lambda_iDlambda_j - lambda_jDlambda_i;
% curl phi = 2*Dlambda_i mycross Dlambda_j;
% psi = lambda_iDlambda_j + lambda_jDlambda_i;
% curl psi = 0;
[Dlambda,volume] = gradbasis3(node,elem);
curlPhi(:,:,6) = 2*mycross(Dlambda(:,:,3),Dlambda(:,:,4),2);
curlPhi(:,:,1) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,2),2);
curlPhi(:,:,2) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,3),2);
curlPhi(:,:,3) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,4),2);
curlPhi(:,:,4) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,3),2);
curlPhi(:,:,5) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,4),2);

%% Assemble matrices
ii = zeros(21*NT,1); jj = zeros(21*NT,1); 
sA = zeros(21*NT,1); sMphi = zeros(21*NT,1); sMpsi = zeros(21*NT,1);
index = 0;
for i = 1:6
    for j = i:6
        % local to global index map and its sign
        signij = dofSign(:,i).*dofSign(:,j);
        idx = (signij == -1);
        % curl-curl matrix
        Aij = dot(curlPhi(:,:,i),curlPhi(:,:,j),2).*volume./epsilon;
        Aij(idx) = -Aij(idx);
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));
        sA(index+1:index+NT) = Aij;
        % mass matrix
        i1 = locEdge(i,1); i2 = locEdge(i,2);
        j1 = locEdge(j,1); j2 = locEdge(j,2);
        % (phi_i,phi_j)
        Mij = 1/20*volume.*( ...
              (1+(i1==j1))*dot(Dlambda(:,:,i2),Dlambda(:,:,j2),2)...
            - (1+(i1==j2))*dot(Dlambda(:,:,i2),Dlambda(:,:,j1),2)...
            - (1+(i2==j1))*dot(Dlambda(:,:,i1),Dlambda(:,:,j2),2)...
            + (1+(i2==j2))*dot(Dlambda(:,:,i1),Dlambda(:,:,j1),2));
        Mij = Mij.*mu;
        Mij(idx) = -Mij(idx);
        sMphi(index+1:index+NT) = Mij;
        % (psi_i,psi_j)
        Mij = 1/20*volume.*( ...
              (1+(i1==j1))*dot(Dlambda(:,:,i2),Dlambda(:,:,j2),2)...
            + (1+(i1==j2))*dot(Dlambda(:,:,i2),Dlambda(:,:,j1),2)...
            + (1+(i2==j1))*dot(Dlambda(:,:,i1),Dlambda(:,:,j2),2)...
            + (1+(i2==j2))*dot(Dlambda(:,:,i1),Dlambda(:,:,j1),2));
        Mij = Mij.*mu;
        sMpsi(index+1:index+NT) = Mij;
        % increment of index
        index = index + NT;
    end
end
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
M = sparse([ii(diagIdx); ii(diagIdx)+NE], [jj(diagIdx); jj(diagIdx)+NE],...
           [sMphi(diagIdx); sMpsi(diagIdx)],Ndof,Ndof);
MU = sparse([ii(upperIdx); ii(upperIdx)+NE], [jj(upperIdx); jj(upperIdx)+NE],...
           [sMphi(upperIdx); sMpsi(upperIdx)],Ndof,Ndof);
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
        Mij = 1/20*volume.*( ...
              (1+(i1==j1))*dot(Dlambda(:,:,i2),Dlambda(:,:,j2),2)...
            - (1+(i1==j2))*dot(Dlambda(:,:,i2),Dlambda(:,:,j1),2)...
            + (1+(i2==j1))*dot(Dlambda(:,:,i1),Dlambda(:,:,j2),2)...
            - (1+(i2==j2))*dot(Dlambda(:,:,i1),Dlambda(:,:,j1),2));
        Mij = Mij.*mu.*double(dofSign(:,j));
        ii(index+1:index+NT) = double(elem2dof(:,i))+NE; 
        jj(index+1:index+NT) = double(elem2dof(:,j));
        ss(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
ML = sparse(ii,jj,ss,Ndof,Ndof);
M  = M  + ML + ML';
clear ML
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
    [lambda,w] = quadpts3(option.fquadorder); % quadrature order is 3
    nQuad = size(lambda,1);
    bt = zeros(NT,12);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:) ... 
             + lambda(p,4)*node(elem(:,4),:);
        Jp = pde.J(pxyz)./epsilon;
    %   locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        for k = 1:6
            i = locEdge(k,1); j = locEdge(k,2);
            rhs = double(dofSign(:,k)).*dot(curlPhi(:,:,k),Jp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
            % psi_k = lambda_iDlambda_j + lambda_jDlambda_i;
            % curl psi = 0
        end
    end
    bt = bt.*repmat(volume,1,12);
    f = accumarray([elem2dof(:); elem2dof(:)+NE],bt(:),[Ndof 1]);
end
clear pxy Jp bt rhs phi_k psi_k curlPhi

%% Set up solver
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 1e4  % Direct solver for small size systems
        option.solver = 'direct';
    else            % Multigrid-type  solver for large size systems
        option.solver = 'cg';
    end
end
solver = option.solver;

%% Assembeling corresponding matrices for HX preconditioner
if ~strcmp(solver,'direct') && ~strcmp(solver,'nosolve')
    AP = sparse(N,N);  % AP = - div(epsilon^{-1}grad) + |mu| I
    BP = sparse(N,N);  % BP = - div(|mu|grad)
    for i = 1:4
        for j = i:4
            temp = dot(Dlambda(:,:,i),Dlambda(:,:,j),2).*volume;
            Aij = temp./epsilon;
            Bij = abs(mu).*temp;
            Mij = 1/20*abs(mu).*volume;
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
% isBdEdge = [];
if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N)
    % Dirichlet boundary condition only
    bdFlag = setboundary3(node,elem,'Dirichlet');
end
if ~isempty(bdFlag)
    % Find boundary edges and nodes
    isBdEdge = false(NE,1);
    isBdEdge(elem2dof(bdFlag(:,1) == 1,[4,5,6])) = true;
    isBdEdge(elem2dof(bdFlag(:,2) == 1,[2,3,6])) = true;
    isBdEdge(elem2dof(bdFlag(:,3) == 1,[1,3,5])) = true;
    isBdEdge(elem2dof(bdFlag(:,4) == 1,[1,2,4])) = true;
    bdDof = [find(isBdEdge); NE + find(isBdEdge)];
    bdEdge = edge(isBdEdge,:);
    isBdNode(bdEdge) = true;
end
if any(isBdEdge)
    % modify the matrix to include Dirichlet boundary condition
    bdidx = zeros(Ndof,1); 
    bdidx(bdDof) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    bigAD = T*(A-M)*T + Tbd;
    if ~strcmp(solver,'direct') && ~strcmp(solver,'nosolve')
        % modify the corresponding Poisson matrix
        bdidx = zeros(N,1); 
        bdidx(isBdNode) = 1;
        Tbd = spdiags(bdidx,0,N,N);
        T = spdiags(1-bdidx,0,N,N);
        AP = T*AP*T + Tbd;
        BP = T*BP*T + Tbd;
    end
else
    bigAD = A - M;
    if ~strcmp(solver,'direct') && ~strcmp(solver,'nosolve')
       BP = BP + 1e-8*speye(N);        
    end    
end

%% Part 2: Find boundary edges and modify the load b
g = zeros(Ndof,1);
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
    face = [1 4 3]; face2locdof = [6 2 3];
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
    face = [1 3 2]; face2locdof = [4 1 2];
    if ~isempty(isBdElem)
        bdb = bdfaceintegral1(isBdElem,face,face2locdof);
        g = g + bdb;
    end
    f = f - g;
end
% nonzero Dirichlet boundary condition
u = zeros(Ndof,1);
if ~isempty(bdEdge) && ~isempty(pde.g_D) && ...
   ~(isnumeric(pde.g_D) && all(pde.g_D == 0))
    % else no bddof or g_D = 0 (no modification needed)
    if (isnumeric(pde.g_D) && length(pde.g_D) == Ndof)
        u(bdDof) = pde.g_D(bdDof);
    else
        u(bdDof) = edgeinterpolate1(pde.g_D,node,bdEdge);
    end
    f = f - (A-M)*u;
    f(bdDof) = u(bdDof);
end
%% Remark
% The order of assign Neumann and Dirichlet boundary condition is
% important to get the right setting of the intersection of Dirichlet and
% Neumann faces.

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 1
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if strcmp(solver,'direct')
    % exact solver
    freeDof = [find(~isBdEdge); NE + find(~isBdEdge)];
    tic
    u(freeDof) = bigAD(freeDof,freeDof)\f(freeDof);    
    time = toc; itStep = 0; flag = 2; err = norm(f - bigAD*u);    
    if (nargin > 4) && isfield(option,'printlevel') && (option.printlevel >= 1)
        fprintf('#dof: %8.0u, Direct solver %4.2g \n',length(f),time);
    end
elseif strcmp(solver,'nosolve')
    eqn = struct('A',A,'M',M,'f',f,'g',g,'bigA',bigAD); info = [];
    return;
else
    % option.x0 = edgeinterpolate1(pde.g_D(node),node,edge);
    option.x0 = u;
    [u,flag,itStep,err,time] = mgMaxwell(bigAD,f,AP,BP,node,elem,edge,HB,isBdEdge,option);
end

%% Output
eqn = struct('A',A,'M',M,'f',f,'g',g,'bigA',bigAD);
info = struct('assembleTime',assembleTime,'solverTime',time,'itStep',itStep,...
              'error',err,'flag',flag);

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
    face2locEdge = [2 3; 3 1; 1 2]; % edge of the face [1 2 3]

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
            pidx = face(face2locEdge(s,1))< face(face2locEdge(s,2));
            % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
            % lambda_i is associated to the local index of the face [1 2 3]
            % Dlambda_j is associtated to the index of tetrahedron
            % - when the direction of the local edge s is consistent with the
            %   global oritentation given in the triangulation, 
            %           s(1) -- k(1),  s(2) -- k(2)
            % - otherwise 
            %           s(2) -- k(1),  s(1) -- k(2)
            if pidx
                phi_k = lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                      - lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,1));
            else
                phi_k = lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                      - lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,1));                   
            end
            rhs = double(dofSign(isBdElem,kk)).*dot(phi_k,gNp,2);
            bt(:,s) = bt(:,s) + w(pp)*rhs; % area is included in normal
            idx(:,s) = elem2dof(isBdElem,kk);
            % psi_k = lambda_iDlambda_j + lambda_jDlambda_i;
            if pidx
                psi_k = lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                      + lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,1));
            else
                psi_k = lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                      + lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,1));                   
            end
            rhs = dot(psi_k,gNp,2);
            bt(:,s+3) = bt(:,s+3) + w(pp)*rhs; % area is included in normal
            idx(:,s+3) = elem2dof(isBdElem,kk) + NE;
        end
    end
    %% Distribute to DOF
    bdb = accumarray(idx(:),bt(:),[Ndof 1]);        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end