function [u,edge,eqn,info] = Maxwellsaddle(node,elem,bdFlag,pde,option,varargin)
%% MAXWELLSADDLE Maxwell equation in mixed form: lowest order edge element.
%
% u = Maxwellsaddle(node,elem,bdFlag,pde) produces the lowest order edge
%   element approximation of the electric field of the magnetostatics
%
%                     curl(mu^(-1)curl u)  = J    in \Omega,  
%                                  div  u  = 0
%                                   n x u = n x g_D  on \Gamma_D,
%                     n x (mu^(-1)curl u) = n x g_N  on \Gamma_N.
% 
% based on the weak formulation
%
% (mu^{-1}curl u, curl v) +(grad p ,v)= (J,v) - <n x g_N,v>_{\Gamma_N}.
% (u,grad q)                          = 0                    
% Assume P\in H_0^1.
% The data of the equation is enclosed in the pde structure:
%   - pde.mu      : permeability, i.e., magnetic constant/tensor
%   - pde.J       : current density
%   - pde.g_D     : Dirichlet boundary condition
%   - pde.g_N     : Neumann boundary condition
%
% The mesh is given by (node,elem) and HB is needed for fast solvers. The
% boundary faces is specified by bdFlag; see <a href="matlab:ifem bddoc">bddoc</a>.
%
% The function Maxwell assembes the matrix equation (A-M)*u = b and solves
% it by the direct solver (small size dof <= 2e3) or the HX preconditioned
% Krylov iterative methods (large size dof > 2e3).
% 
% u = Maxwellsaddle(node,elem,bdFlag,pde,option) specifies the solver options.
%   - option.solver == 'direct': the built in direct solver \ (mldivide)
%   - option.solver == 'mg':     multigrid-type solvers mg is used.
%   - option.solver == 'none': the solution u = u_D. 
% The default setting is to use the direct solver for small size problems
% and multigrid solvers for large size problems. For more options on the
% multigrid solver mg, type help mg.
%
% [u,edge] = Maxwellsaddle(node,elem,pde,bdEdge) returns also the edge array
% which is essential for edge elements. 
%
% [u,edge,eqn] = Maxwellsaddle(node,elem,pde,bdEdge) returns also the equation
% structure eqn, which includes: 
% - eqn.A: matrix for differential operator;
% - eqn.M: mass matrix;
% - eqn.f: right hand side 
% - eqn.g: vector enclosed the Neumann boundary condition
%
% [u,edge,eqn,info] = Maxwellsaddle(node,elem,pde,bdEdge) returns also the
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
%   Maxwellsaddlemgtest
%
% See also Maxwell1, Maxwell2, cubeMaxwell, mgMaxwell
%
% Reference page in Help browser
%       <a href="matlab:ifem Maxwelldoc">Maxwelldoc</a> 
%
% Created by Jie Zhou based on Maxwell(node,elem,pde,bdEdge) on
% 09,Sep,2013.
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.


%% Set up optional input arguments
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end

t = cputime;
%% Sort elem to ascend ordering
elemMG  = elem;     % save elem and bdFlag for multigrid
bdFlagMG = bdFlag;
[elem,bdFlag] = sortelem3(elem,bdFlag);
[elem2dof,edge] = dof3edge(elem);
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
N = size(node,1);   NT = size(elem,1);  Ne = size(edge,1);

%% Compute coefficients
if ~isfield(pde,'mu'), pde.mu = 1; end
if ~isempty(pde.mu) && isnumeric(pde.mu)
    mu = pde.mu;                % mu is an array
else                            % mu is a function
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;
    mu = pde.mu(center);              
end
if ~isfield(pde,'epsilon'), pde.epsilon = 0; end
if ~isempty(pde.epsilon) && isnumeric(pde.epsilon)
    epsilon = pde.epsilon;      % epsilon is an array
else                            % epsilon is a function
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
              node(elem(:,3),:) + node(elem(:,4),:))/4;
    epsilon = pde.epsilon(center);              
end
if isfield(pde,'omega')
    omega = pde.omega;
else
    omega = 1;
end
epsilon = omega^2*epsilon; 

%% Element-wise basis
% edge indices of 6 local bases: 
% [1 2], [1 3], [1 4], [2 3], [2 4], [3 4]
% phi = lambda_i Dlambda_j - lambda_j Dlambda_i;
% curl phi = 2*Dlambda_i \times Dlambda_j;
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
sA = zeros(21*NT,1); sM = zeros(21*NT,1);
index = 0;
for i = 1:6
    for j = i:6
        % local to global index map
        % curl-curl matrix
        Aij = dot(curlPhi(:,:,i),curlPhi(:,:,j),2).*volume./mu;
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));
        sA(index+1:index+NT) = Aij;
        % mass matrix
        % locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        i1 = locEdge(i,1); i2 = locEdge(i,2);
        j1 = locEdge(j,1); j2 = locEdge(j,2);
        Mij = 1/20*volume.*( (1+(i1==j1))*DiDj(:,i2,j2) ...
                           - (1+(i1==j2))*DiDj(:,i2,j1) ...
                           - (1+(i2==j1))*DiDj(:,i1,j2) ...
                           + (1+(i2==j2))*DiDj(:,i1,j1));
        %Mij = Mij.*epsilon;
        sM(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
clear curlPhi % clear large size data
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ne,Ne);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ne,Ne);
A = A + AU + AU';
M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ne,Ne);
MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ne,Ne);
M = M + MU + MU';
grad =icdmat(double(edge),[-1,1]);
% G    = grad'*M;
G    = M*grad;

%% Assemble right hand side
f = zeros(Ne,1);
if ~isfield(pde,'J') || (isfield(pde,'J') && isreal(pde.J) && all(pde.J==0))
    pde.J = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order is 3
end
if isfield(pde,'J') && ~isempty(pde.J) && ~isnumeric(pde.J)
    [lambda,w] = quadpts3(option.fquadorder);
    nQuad = size(lambda,1);
    bt = zeros(NT,6);
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
        end
    end
    bt = bt.*repmat(volume,1,6);
    f = accumarray(elem2dof(:),bt(:),[Ne 1]);
elseif isfield(pde,'J') && ~isempty(pde.J) && isnumeric(pde.J)
    switch size(pde.J,1)
        case Ne % rhs already computed
            f = pde.J;
        case NT % piecwise constant
            bt = zeros(NT,6);
            for k = 1:6
                i = locEdge(k,1); j = locEdge(k,2);
                % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
                phi_k = (Dlambda(:,:,j)-Dlambda(:,:,i))/4;
                bt(:,k) = dot(phi_k,pde.J,2);
            end
            bt = bt.*repmat(volume,1,6);
            f = accumarray(elem2dof(:),bt(:),[Ne 1]);
    end
end
clear pxyz Jp bt rhs phi_k

%% Record assembling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 1
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end
%% Boundary condtiions
[AD,f,u,freeNode,freeEdge] = getbdMaxwellsaddle(f);
p  = zeros(N,1);
g0 = -G'*u;
% reduce to free dofs
A  = A(freeEdge,freeEdge);
G  = G(freeEdge,freeNode);
grad = grad(freeEdge,freeNode);
M = M(freeEdge,freeEdge);
f = f(freeEdge);
g0 = g0(freeNode);

%% Solve the linear system
if isfield(option,'solver')
    method = upper(option.solver);
else
    method = 'TRI';
end

% multigrid options 
option.mg.isFreeEdge = freeEdge; % needed in mg solver
option.mg.isFreeNode = freeNode; 
switch method
    case 'DIRECT'
        t = cputime;
        Ni = size(freeNode,1); Nei = size(freeEdge,1);
        bigA   = [A  G; ...
                  G' sparse(Ni,Ni)];
        temp   = bigA\[f; g0];
        u(freeEdge) = temp(1:Nei);
        p(freeNode) = temp(Nei+1:Nei+Ni);
        residual = norm([f;g0] - bigA*temp);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'MG'
        [u0,p0,info] = mgMaxwellsaddle(A,G,f,g0,node,elemMG,bdFlagMG,M,grad,option.mg);
        u(freeEdge)  = u0;
        p(freeNode)  = p0;        
    case 'TRI' % GMRES with a triangular preconditioner 
        [u0,p0,info] = tripreMaxwellsaddle(A,G,f,g0,node,elemMG,bdFlagMG,M,grad,option.mg);
        u(freeEdge)  = u0;
        p(freeNode)  = p0;
    case 'DIAG' % MINRES with a diagonal preconditioner
        [u0,p0,info] = diapreMaxwellsaddle(A,G,f,g0,node,elemMG,bdFlagMG,option.mg);
        u(freeEdge)  = u0;
        p(freeNode)  = p0;
end

%check the Langrange is correct or not.
% normp  = norm(p);
% if(normp>1.0/N)
% disp('the Langrange multer is wrong')
% end

%% Output
eqn = struct('A',A,'M',M,'f',f,'g',g0,'bigA',AD,'freeEdge',freeEdge);
info.assembleTime = assembleTime;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdMaxwellsaddle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,f,u,freeNode,freeEdge,isPureNeumann] = getbdMaxwellsaddle(f)

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
        isBdEdge = false(Ne,1);
        isBdNode = false(N,1);
        isBdEdge(elem2dof(bdFlag(:,1) == 1,[4,5,6])) = true;
        isBdEdge(elem2dof(bdFlag(:,2) == 1,[2,3,6])) = true;
        isBdEdge(elem2dof(bdFlag(:,3) == 1,[1,3,5])) = true;
        isBdEdge(elem2dof(bdFlag(:,4) == 1,[1,2,4])) = true;
        bdEdge = edge(isBdEdge,:);
        isBdNode(bdEdge) = true;
    end
    % modify the matrix to include the Dirichlet boundary condition
    if any(isBdEdge)  % contains Dirichlet boundary condition
        bdidx = zeros(Ne,1); 
        bdidx(isBdEdge) = 1;
        Tbd = spdiags(bdidx,0,Ne,Ne);
        Te = spdiags(1-bdidx,0,Ne,Ne);
        AD = Te*(A - epsilon*M)*Te + Tbd;
    else      % pure Neumann boundary condition
        isPureNeumann = 1;
        AD = A - epsilon*M  ;
    end

    %% Part 2: Find boundary edges and modify the load b
    g = zeros(Ne,1);
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
            bdb = bdfaceintegral(isBdElem,face,face2locdof);
            g = bdb;
        end
        % face 2
        isBdElem = find(bdFlag(:,2) == 2);
        face = [1 4 3]; face2locdof = [6 2 3];
        if ~isempty(isBdElem)
            bdb = bdfaceintegral(isBdElem,face,face2locdof);
            g = g + bdb; 
        end
        % face 3
        isBdElem = find(bdFlag(:,3) == 2);
        face = [1 2 4]; face2locdof = [5 3 1];
        if ~isempty(isBdElem)
            bdb = bdfaceintegral(isBdElem,face,face2locdof);
            g = g + bdb; 
        end
        % face 4
        isBdElem = find(bdFlag(:,4) == 2);
        face = [1 3 2]; face2locdof = [4 1 2];
        if ~isempty(isBdElem)
            bdb = bdfaceintegral(isBdElem,face,face2locdof);
            g = g + bdb;
        end
        f = f - g;
    end
    % nonzero Dirichlet boundary condition
    u = zeros(Ne,1);
    if ~isempty(bdEdge) && ~isempty(pde.g_D) && ...
       ~(isnumeric(pde.g_D) && all(pde.g_D == 0))
        % else no bddof or g_D = 0 (no modification needed)
        if (isnumeric(pde.g_D) && length(pde.g_D) == Ne)
            u(isBdEdge) = pde.g_D(isBdEdge);
        else
            u(isBdEdge) = edgeinterpolate(pde.g_D,node,bdEdge);
        end
        f = f - (A - epsilon*M)*u;
        f(isBdEdge) = u(isBdEdge);
    end
    % find free node, free edge and free dof
    freeEdge = find(~isBdEdge);
    freeNode = find(~isBdNode);
    %% Remark
    % The order of assign Neumann and Dirichlet boundary condition is
    % important to get the right setting of the intersection of Dirichlet and
    % Neumann faces.
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions bdfaceintegral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function bdb = bdfaceintegral(isBdElem,face,face2locdof)
    %% Compute boundary surface integral of lowest order edge element.
    %  bdb(k) = \int_{face} (nï¿?g_N, phi_k) dS

    %% Compute scaled normal
    faceIdx = true(4,1);
    faceIdx(face) = false;
    normal = -3*repmat(volume(isBdElem),1,3).*Dlambda(isBdElem,:,faceIdx);

    %% Data structure
    tetLocEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % edge of a tetrahedral [1 2 3 4]
    face2locEdge = [2 3; 3 1; 1 2]; % edge of the face [1 2 3]

    %% Compute surface integral
    Nbd = length(isBdElem);
    bt = zeros(Nbd,3);
    idx = zeros(Nbd,3,'int32');
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
            rhs = dot(phi_k,gNp,2);
            bt(:,s) = bt(:,s) + w(pp)*rhs; % area is included in normal; see line 28
            idx(:,s) = elem2dof(isBdElem,kk);
        end
    end
    %% Distribute to DOF
    bdb = accumarray(idx(:),bt(:),[Ne 1]);        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end