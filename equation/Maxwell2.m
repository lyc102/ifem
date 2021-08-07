function [u,T,eqn,info] = Maxwell2(node,elem,bdFlag,pde,option,varargin)
%% MAXWELL2 Maxwell equation: quadratic edge element.
%
% u = Maxwell2(node,elem,HB,pde,bdFlag) produces the quadratic (1st kind)
% edge element approximation of the electric field of the time harmonic
% Maxwell equation.
%
% curl(mu^(-1)curl u) + omega*epsilon u = J    in \Omega,  
%                                 n × u = n × g_D  on \Gamma_D,
%                   n × (mu^(-1)curl u) = n × g_N  on \Gamma_N.
% 
% based on the weak formulation
%
% (mu^{-1}curl u, curl v) + (epsilon u,v) = (J,v) - <n × g_N,v>_{\Gamma_N}.
%
% The data of the equation is enclosed in the pde structure:
%   - pde.mu      : permeability, i.e., magnetic constant/tensor
%   - pde.epsilon : a complex dielectric constant/tensor
%   - pde.omega   : wave number
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
% u = Maxwell(node,elem,HB,pde,bdFlag,option) specifies the solver options.
%   - option.solver == 'direct': the built in direct solver \ (mldivide)
%   - option.solver == 'mg':     multigrid-type solvers mg is used.
%   - option.solver == 'notsolve': the solution u = u_D. 
% The default setting is to use the direct solver for small size problems
% and multigrid solvers for large size problems. For more options on the
% multigrid solver mg, type help mg.
%
% [u,T] = Poisson(node,elem,pde,bdEdge) returns also the structure T containing 
% edge, face information of the current triangulation
%
% [u,T,eqn] = Poisson(node,elem,pde,bdEdge) returns also the equation
% structure eqn, which includes: 
% - eqn.A: matrix for differential operator;
% - eqn.M: mass matrix;
% - eqn.f: right hand side 
% - eqn.g: vector enclosed the Neumann boundary condition
%
%
% Example
%
%   cubeMaxwell2
%
% See also Maxwell, Maxwell1, cubeMaxwell2, mgMaxwell2
%
% <a href="matlab:ifem Maxwell2doc">Maxwell2 doc</a> 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

tstart = cputime;
%% Set up optional input arguments
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
if nargin>=6
    HB = varargin{1};
else
    HB = [];
end        

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
    epsilon = pde.omega*epsilon; 
end

%% Sort elem to ascend ordering
elemold = elem;
[elem,bdFlag] = sortelem3(elem,bdFlag);

%% Construct Data Structure
[Dlambda,volume] = gradbasis3(node,elem);
% locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
% locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
locBasesIdx = [1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % phi
               1 2 0; 1 3 0; 1 4 0; 2 3 0; 2 4 0; 3 4 0; ... % psi
               3 2 4; 3 1 4; 2 1 4; 2 1 3; ...
               4 2 3; 4 1 3; 4 1 2; 3 1 2]; % chi
% See <matlab:ifem('Maxwell2doc') Maxwell2doc> for details.            
% construct elem2edge and elem2face 
[elem2edge,edge] = dof3edge(elem);
[elem2face,face] = dof3face(elem); 
N = size(node,1);  NT = size(elem,1);  
NE = size(edge,1); NF = size(face,1); Ndof = 2*(NE + NF);
elem2dof = [elem2edge elem2edge+NE elem2face+2*NE elem2face+2*NE+NF];
face2edge = zeros(NF,3,'int32');
% face is given by auxtructure3 and is sorted according to global indices.
% locFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
% face2edge is used to compute uI. So it is consistent with the local index
% system in edgeinterpolate2, i.e., if face is (i,j,k) with i<j<k, then the
% three edges are [i j], [i k], [j k]. 
face2edge(elem2face(:,1),:) = elem2edge(:,[4 5 6]);
face2edge(elem2face(:,2),:) = elem2edge(:,[2 3 6]);
face2edge(elem2face(:,3),:) = elem2edge(:,[1 3 5]);
face2edge(elem2face(:,4),:) = elem2edge(:,[1 2 4]);

%% Triangulation output
T = struct('edge',edge,'face',face,'face2edge',face2edge, ...
           'elem',elem,'bdFlag',bdFlag);
 
%% Assemble matrices for second family linear element
DiDjcross = zeros(NT,3,4,4);
for i = 1:4
    for j = i+1:4        
        DiDjcross(:,:,i,j) = mycross(Dlambda(:,:,i),Dlambda(:,:,j),2);
        DiDjcross(:,:,j,i) = -DiDjcross(:,:,i,j);
    end
end
DiDj = zeros(NT,4,4);
for i = 1:4
    for j = i:4        
        DiDj(:,i,j) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2);
        DiDj(:,j,i) = DiDj(:,i,j);
    end
end

ii = zeros(210*NT,1); jj = zeros(210*NT,1); 
index = 0;
[lambda,w] = quadpts3(3); % quadrature order is 3
nQuad = size(lambda,1);
for i = 1:20
    for j = i:20
        ii(index+1:index+NT) = double(elem2dof(:,i));
        jj(index+1:index+NT) = double(elem2dof(:,j));
        i1 = locBasesIdx(i,1); i2 = locBasesIdx(i,2); i3 = locBasesIdx(i,3);
        j1 = locBasesIdx(j,1); j2 = locBasesIdx(j,2); j3 = locBasesIdx(j,3);
        Aij = zeros(NT,1);  Mij = zeros(NT,1); 
        %% curl-curl matrix
        if (i<=6) && (j<=6)
            Aij = 4*dot(DiDjcross(:,:,i1,i2),DiDjcross(:,:,j1,j2),2);
        end
        if (i<=6) && (j>12)
            Aij = dot(DiDjcross(:,:,i1,i2),...
            DiDjcross(:,:,j1,j3)-DiDjcross(:,:,j1,j2)+2*DiDjcross(:,:,j2,j3),2)/2;
        end
        if (i>12) && (j>12)
        % curl chi_i =  Dlambda_{i1}mycross phi_{i2,i3} + lambda_{i1}curl phi_{i2,i3}
        % curl chi_j =  Dlambda_{j1}mycross phi_{j2,j3} + lambda_{j1}curl phi_{j2,j3}
            % (Dlambda_{i1}mycross phi_{i2,i3}) dot (Dlambda_{j1}mycross phi_{j2,j3})
            temp11 = ((1+(i2==j2))*dot(DiDjcross(:,:,i1,i3),DiDjcross(:,:,j1,j3),2) ...
                    - (1+(i2==j3))*dot(DiDjcross(:,:,i1,i3),DiDjcross(:,:,j1,j2),2) ...
                    - (1+(i3==j2))*dot(DiDjcross(:,:,i1,i2),DiDjcross(:,:,j1,j3),2) ...
                    + (1+(i3==j3))*dot(DiDjcross(:,:,i1,i2),DiDjcross(:,:,j1,j2),2))/20;
            % lambda_{i1}curl phi_{i2,i3} dot lambda_{j1}curl phi_{j2,j3}
            temp22 = (1+(i1==j1))/5*dot(DiDjcross(:,:,i2,i3),DiDjcross(:,:,j2,j3),2);
            % Dlambda_{i1}mycross phi_{i2,i3} dot lambda_{j1}curl phi_{j2,j3}
            temp12 = dot( (1+(j1==i2))*DiDjcross(:,:,i1,i3) ...
                        - (1+(j1==i3))*DiDjcross(:,:,i1,i2), ...
                                       DiDjcross(:,:,j2,j3),2)/10;                         
            % Dlambda_{j1}mycross phi_{j2,j3} dot lambda_{i1}curl phi_{i2,i3}
            temp21 = dot( (1+(i1==j2))*DiDjcross(:,:,j1,j3) ...
                        - (1+(i1==j3))*DiDjcross(:,:,j1,j2), ...
                                       DiDjcross(:,:,i2,i3),2)/10;   
            Aij = temp11 + temp22 + temp12 + temp21;
        end
        if (i<=6) && (j<=6)
            % block 1: (phi_i,phi_j)
            Mij = 1/20*((1+(i1==j1))*DiDj(:,i2,j2) ...
                      - (1+(i1==j2))*DiDj(:,i2,j1) ...
                      - (1+(i2==j1))*DiDj(:,i1,j2) ...
                      + (1+(i2==j2))*DiDj(:,i1,j1));
        end
        if (i<=6) && (7<=j) && (j<=12)
            % block 2: (psi_j,phi_i)
            Mij = 1/20*( (1+(j1==i1))*DiDj(:,j2,i2) ...
                       - (1+(j1==i2))*DiDj(:,j2,i1) ...
                       + (1+(j2==i1))*DiDj(:,j1,i2) ...
                       - (1+(j2==i2))*DiDj(:,j1,i1));

        end
        if (7<=i) && (i<=12) && (7<=j) && (j<=12)
            % block 3: (psi_j,psi_i)
            Mij = 1/20*((1+(i1==j1))*DiDj(:,i2,j2) ...
                      + (1+(i1==j2))*DiDj(:,i2,j1) ...
                      + (1+(i2==j1))*DiDj(:,i1,j2) ...
                      + (1+(i2==j2))*DiDj(:,i1,j1));
        end
        if (i<=6) && (j>12)
            % block 4: (chi_j,phi_i)
            Mij = intlambda([j1,i1,j2],3)*DiDj(:,i2,j3) ...
                 -intlambda([j1,i1,j3],3)*DiDj(:,i2,j2) ...
                 -intlambda([j1,i2,j2],3)*DiDj(:,i1,j3) ...
                 +intlambda([j1,i2,j3],3)*DiDj(:,i1,j2);            
        end
        if (7<=i) && (i<=12) && (j>12)
            % block 5: (chi_j,psi_i)
            Mij = intlambda([j1,i1,j2],3)*DiDj(:,i2,j3) ...
                 -intlambda([j1,i1,j3],3)*DiDj(:,i2,j2) ...
                 +intlambda([j1,i2,j2],3)*DiDj(:,i1,j3) ...
                 -intlambda([j1,i2,j3],3)*DiDj(:,i1,j2);                        
        end
        if (i>12) && (j>12)
            % block 6: (chi_j,chi_i)
            Mij = intlambda([i1,j1,i2,j2],3)*DiDj(:,i3,j3) ...
                 -intlambda([i1,j1,i2,j3],3)*DiDj(:,i3,j2) ...
                 -intlambda([i1,j1,i3,j2],3)*DiDj(:,i2,j3) ...
                 +intlambda([i1,j1,i3,j3],3)*DiDj(:,i2,j2);                        
        end
        Aij = Aij.*volume./mu;
        Mij = Mij.*volume.*epsilon;
        sA(index+1:index+NT) = Aij;
        sM(index+1:index+NT) = Mij;
        index = index + NT;
    end
end
clear curlBasis_i curlBasis_j basis_i basis_j
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
M = sparse(ii(diagIdx),jj(diagIdx),sM(diagIdx),Ndof,Ndof);
MU = sparse(ii(upperIdx),jj(upperIdx),sM(upperIdx),Ndof,Ndof);
M = M + MU + MU';
clear AU MU

%% Assemble right hand side
f = zeros(Ndof,1);
if ~isfield(pde,'J') || (isfield(pde,'J') && isreal(pde.J) && all(pde.J==0))
    pde.J = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 4;   % default order is 3
end
if isfield(pde,'J') && ~isempty(pde.J)
    [lambda,w] = quadpts3(option.fquadorder); % quadrature order is 4
    nQuad = size(lambda,1);
    bt = zeros(NT,20);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:) ... 
             + lambda(p,4)*node(elem(:,4),:);
        Jp = pde.J(pxyz);
        for k = 1:20
            k1 = locBasesIdx(k,1); 
            k2 = locBasesIdx(k,2); 
            k3 = locBasesIdx(k,3);
            % evaluate basis at quadrature points
            if k<=6
            % phi_k = lambda_{k1}Dlambda_{k2} - lambda_{k2}Dlambda_{k1};
                basis_k = (lambda(p,k1)*Dlambda(:,:,k2) ...
                          -lambda(p,k2)*Dlambda(:,:,k1));
            elseif k<=12
            % phi_k = lambda_{k1}Dlambda_{k2} + lambda_{k2}Dlambda_{k1};
                basis_k = (lambda(p,k1)*Dlambda(:,:,k2) ...
                          +lambda(p,k2)*Dlambda(:,:,k1));
            else
            % chi_k = lambda_{k1}phi_{k2,k3};    
                basis_k = lambda(p,k1)*(lambda(p,k2)*Dlambda(:,:,k3) ...
                                       -lambda(p,k3)*Dlambda(:,:,k2));
            end            
            rhs = dot(basis_k,Jp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end
    bt = bt.*repmat(volume,1,20);
    f = accumarray(elem2dof(:),bt(:),[Ndof 1]);
end
clear pxy Jp bt rhs basis_k


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
if ~strcmp(solver,'direct')
    AP = sparse(N,N);  % AP = - div(mu^{-1}grad) + |epsilon| I
    BP = sparse(N,N);  % BP = - div(|epsilon|grad)
    for i = 1:4
        for j = i:4
            temp = dot(Dlambda(:,:,i),Dlambda(:,:,j),2).*volume;
            Aij = 1./mu.*temp;
            Bij = abs(real(epsilon)).*temp;
            Mij = 1/20*abs(real(epsilon)).*volume;
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
if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N)
    % Dirichlet boundary condition only
    bdFlag = setboundary3(node,elem,'Dirichlet');
end
isBdDof = false(Ndof,1);
if ~isempty(bdFlag)
    %% Dirichlet boundary condition on edge dofs
    % Find boundary faces, edges and nodes
    isBdFace = false(NF,1);
    isBdFace(elem2face(bdFlag(:,1) == 1,1)) = true;
    isBdFace(elem2face(bdFlag(:,2) == 1,2)) = true;
    isBdFace(elem2face(bdFlag(:,3) == 1,3)) = true;
    isBdFace(elem2face(bdFlag(:,4) == 1,4)) = true;
    boundaryFace = face(isBdFace,:);
    bdFace2edge = face2edge(isBdFace,:);
    isBdEdge = false(NE,1);
    isBdEdge(bdFace2edge(:)) = true;
    edgeBdDof = [find(isBdEdge); NE + find(isBdEdge)];
    bdEdge = edge(isBdEdge,:);
    isBdNode(bdEdge) = true;
%     bdNode = find(isBdNode);
    faceBdDof = 2*NE + [find(isBdFace); NF+find(isBdFace)];
    edgeIdxMap = zeros(NE,1);
    edgeIdxMap(isBdEdge) = 1:size(bdEdge,1);
    bdFace2edge = edgeIdxMap(bdFace2edge);
    isBdDof(edgeBdDof) = true;
    isBdDof(faceBdDof) = true;
end
if any(isBdEdge)
    % modify the matrix to include Dirichlet boundary condition
    bdidx = zeros(Ndof,1); 
    bdidx(isBdDof) = 1;
    Pbd = spdiags(bdidx,0,Ndof,Ndof);
    P = spdiags(1-bdidx,0,Ndof,Ndof);
    bigAD = P*(A-M)*P + Pbd;
    if ~strcmp(solver,'direct')
        % modify the corresponding Poisson matrix
        bdidx = zeros(N,1); 
        bdidx(isBdNode) = 1;
        Pbd = spdiags(bdidx,0,N,N);
        P = spdiags(1-bdidx,0,N,N);
        AP = P*AP*P + Pbd;
        BP = P*BP*P + Pbd;
    end
else
    bigAD = A - M;
    if ~strcmp(solver,'direct')
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
        bdb = bdfaceintegral2(isBdElem,face,face2locdof);
        g = bdb; 
    end
    % face 2
    isBdElem = find(bdFlag(:,2) == 2);
    face = [1 3 4]; face2locdof = [6 3 2];
    if ~isempty(isBdElem)
        bdb = bdfaceintegral2(isBdElem,face,face2locdof);
        g = g + bdb;
    end
    % face 3
    isBdElem = find(bdFlag(:,3) == 2);
    face = [1 2 4]; face2locdof = [5 3 1];
    if ~isempty(isBdElem)
        bdb = bdfaceintegral2(isBdElem,face,face2locdof);
        g = g + bdb;
    end
    % face 4
    isBdElem = find(bdFlag(:,4) == 2);
    face = [1 2 3]; face2locdof = [4 2 1];
    if ~isempty(isBdElem)
        bdb = bdfaceintegral2(isBdElem,face,face2locdof);
        g = g + bdb;
    end
    f = f - g;
end
% nonzero Dirichlet boundary condition
u = zeros(Ndof,1);
if ~isempty(bdEdge) && ~isempty(pde.g_D) && ...
   ~(isnumeric(pde.g_D) && all(pde.g_D == 0))
    % else no bddof or g_D = 0 (no modification needed)
    % compute boundary values of edge dofs
    if (isnumeric(pde.g_D) && length(pde.g_D) > 1)
        pde.g_D(end+1:Ndof) = 0;
    end
    if (isnumeric(pde.g_D) && length(pde.g_D) == Ndof)
        u(isBdDof) = pde.g_D(isBdDof);
    else
        u(isBdDof) = edgeinterpolate2(pde.g_D,node,bdEdge,...
                                  boundaryFace,bdFace2edge);
    end
    f = f - (A-M)*u;
    % modify the matrix to include Dirichlet boundary condition
    f(isBdDof) = u(isBdDof);
end
%% Remark
% The order of assign Neumann and Dirichlet boundary condition is
% important to get the right setting of the intersection of Dirichlet and
% Neumann faces.

%% Record assembling time
assembleTime = cputime - tstart;
if ~isfield(option,'printlevel'), option.printlevel = 1; end

%% Solve the system of linear equations
if strcmp(solver,'direct')
    % exact solver
    freeDof = find(~isBdDof);
    u(freeDof) = bigAD(freeDof,freeDof)\f(freeDof);
    time = cputime - tstart; 
    err = norm(f - bigAD*u);
    info = struct('solverTime',time,'assembleTime',assembleTime,'itStep',0, ...
                  'stopErr',0, 'error',err,'flag',2);      
    if (nargin > 4) && isfield(option,'printlevel') && (option.printlevel >= 1)
        fprintf('#dof: %8.0u, Direct solver %4.2g \n',length(f),time);
    end
elseif strcmp(solver,'none')
    eqn = struct('A',A,'M',M,'AP',AP,'BP',BP,'f',f,'g',g,'bigA',bigAD,'isBdEdge',isBdEdge); 
    info = [];
    return;
elseif strcmp(solver,'amg')
%     u0 = edgeinterpolate(pde.g_D,node,edge);
    u0 = u;
    option.x0 = u0;
    option.alpha = ones(NE,1);
    option.beta = ones(NE,1);
%     option.isBdEdge = isBdEdge;
    option.outsolver = 'cg';    
    [u,info] = amgMaxwell(bigAD,f,node,edge,option);
else
%    option.x0 = edgeinterpolate2(pde.g_D,node,T.edge,T.face,T.face2edge);
    option.x0 = u;
    option.outsolver = 'cg';    
    [u,info] = mgMaxwell(bigAD,f,AP,BP,node,elemold,edge,HB,isBdEdge,option);
end

%% Output
eqn = struct('A',A,'M',M,'f',f,'g',g,'bigA',bigAD);
info.assembleTime = assembleTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions bdfaceintegral2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function bdb = bdfaceintegral2(isBdElem,face,face2locdof)
    %% Compute boundary surface integral of quadratic edge element.
    %  bdb(k) = \int_{face} (n×g_N, phi_k) dS
    %  bdb(k+3) = \int_{face} (n×g_N, psi_k) dS
    %  bdb(k+6) = \int_{face} (n×g_N, chi_k) dS
    
    %% Compute scaled normal
    faceIdx = true(4,1);
    faceIdx(face) = false;
    normal = -3*repmat(volume(isBdElem),1,3).*Dlambda(isBdElem,:,faceIdx);
    
    %% Data structure
    tetLocEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % edge of a tetrahedral [1 2 3 4]
    face2locEdge = [2 3; 1 3; 1 2]; % edge of the face [1 2 3]

    %% Compute surface integral
    Nbd = length(isBdElem);
    bt = zeros(Nbd,8);
    idx = zeros(Nbd,8,'int32');
    [lambda,w] = quadpts(3); % quadrature order is 3
    nQuad = size(lambda,1);
    for pp = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(pp,1)*node(elem(isBdElem,face(1)),:) ...
             + lambda(pp,2)*node(elem(isBdElem,face(2)),:) ... 
             + lambda(pp,3)*node(elem(isBdElem,face(3)),:);
        gNp = pde.g_N(pxyz,normal);    
        % basis associated to edges
        for s = 1:3
            kk = face2locdof(s);
            phi_k = lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                  - lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,1));
            rhs = dot(phi_k,gNp,2);
            bt(:,s) = bt(:,s) + w(pp)*rhs; % area is included in normal; see line 28
            idx(:,s) = elem2dof(isBdElem,kk);
            % psi_k = lambda_iDlambda_j + lambda_jDlambda_i;
            psi_k = lambda(pp,face2locEdge(s,1))*Dlambda(isBdElem,:,tetLocEdge(kk,2)) ...
                  + lambda(pp,face2locEdge(s,2))*Dlambda(isBdElem,:,tetLocEdge(kk,1));
            rhs = dot(psi_k,gNp,2);
            bt(:,s+3) = bt(:,s+3) + w(pp)*rhs; % area is included in normal; see line 27
            idx(:,s+3) = elem2dof(isBdElem,kk) + NE;
        end
        % basis associated to faces
        % chi_1 = lambda_{2}phi_{1,3};
        % chi_2 = lambda_{3}phi_{1,2};
        chi_1 = lambda(pp,2)*(lambda(pp,1)*Dlambda(isBdElem,:,face(3)) ...
                             -lambda(pp,3)*Dlambda(isBdElem,:,face(1)));
        rhs = dot(chi_1,gNp,2);
        bt(:,7) = bt(:,7) + w(pp)*rhs;
        idx(:,7) = elem2dof(isBdElem,12+find(faceIdx));
        chi_2 = lambda(pp,3)*(lambda(pp,1)*Dlambda(isBdElem,:,face(2)) ...
                             -lambda(pp,2)*Dlambda(isBdElem,:,face(1)));
        rhs = dot(chi_2,gNp,2);
        bt(:,8) = bt(:,8) + w(pp)*rhs;
        idx(:,8) = elem2dof(isBdElem,12+find(faceIdx)+4);
    end
    %% Distribute to DOF
    bdb = accumarray(idx(:),bt(:),[Ndof 1]);        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end