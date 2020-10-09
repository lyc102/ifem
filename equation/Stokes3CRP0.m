function [soln,eqn,info] = Stokes3CRP0(node,elem,bdFlag,pde,option,varargin)
%% Stokes3CRP0 solves Stokes equation
%
%   [u,p] = STOKESPCR(node,elem,bdFlag,pde) use Crouzeix and Raviart
%   nonconforming elements to approximate velocity u and piecewise constant
%   to approximate pressure p, repectively.
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0 in \Omega,
%   with
%       Dirichlet boundary condition        u = g_D  on \Gamma_D,
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.
% 
%   Update notes: (April 2020) all the DoF indexing is now boolean/logical.
%  
% See also Stokes, StokesP2P1, PoissonCR.
%
% Add by Shuhao Cao. Apr, 2020. 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end
%% Construct Data Structure
[elem2face,face] = dof3face(elem);
NF = size(face,1); NT = size(elem,1); Nu = NF; Np = NT;

t = cputime;
%% Compute geometric quantities and gradient of local basis
[Dlambda,volume] = gradbasis3(node,elem);

%% Assemble stiffness matrix for Laplace operator
A = sparse(Nu,Nu);
for i = 1:4
    for j = i:4
        % local to global index map
        ii = double(elem2face(:,i));
        jj = double(elem2face(:,j));
        % local stiffness matrix
        Aij = 9*dot(Dlambda(:,:,i),Dlambda(:,:,j),2).*volume;
        if (j==i)
            A = A + sparse(ii,jj,Aij,Nu,Nu);
        else
            A = A + sparse([ii,jj],[jj,ii],[Aij; Aij],Nu,Nu);
        end
    end
end
clear Aij
A = blkdiag(A,A,A);

%% Assemble matrix for divergence operator
Dpsi = -3*reshape(permute(Dlambda, [1 3 2]), [4*NT, 3])...
            .*repmat(volume,[4,3]);

B = sparse(NT, 3*NF);
for j = 1:3 % components of \partial_{x_j}
    idxj = double(elem2face(:)) + (j-1)*NF;
    B = B + sparse(repmat((1:NT)',[4,1]), idxj, Dpsi(:,j),NT,3*NF);
end

%% Assemble right hand side
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end

if ~isfield(pde,'f') || (isreal(pde.f) && all(pde.f==0, 'all'))
    pde.f = [];
end

[lambda,weight] = quadpts3(option.fquadorder);
nQuad = size(lambda,1);
fxyz = zeros(Nu,3); % x,y,z components for RHS
if ~isempty(pde.f) && ~isnumeric(pde.f)
    for j = 1:3 % 3 components, each has NF DoFs
        ftpsi = zeros(NT,4);
        for p = 1:nQuad
            % quadrature points in the x-y-z coordinate
            pxyz = lambda(p,1)*node(elem(:,1),:) ...
                 + lambda(p,2)*node(elem(:,2),:) ...
                 + lambda(p,3)*node(elem(:,3),:) ...
                 + lambda(p,4)*node(elem(:,4),:);
            fp = pde.f(pxyz);            
            for i = 1:4 % 4 faces
                ftpsi(:,i) = ftpsi(:,i) + weight(p)*(1-3*lambda(p,i))*fp(:,j).*volume;
            end
        end
        fxyz(:,j) = accumarray(elem2face(:),ftpsi(:),[Nu 1]);
    end
elseif ~isempty(pde.f) && isnumeric(pde.f) % f is a polynomial vector
    switch size(pde.f,1)
        case NT % piecewise constant vector
            for j = 1:3 % 3 components
                for i = 1:4 % 4 faces
                    fxyz(:,j) = fxyz(:,j) + ...
                        accumarray(elem2face(:,i), pde.f(:,j).*volume/4, [NF 1]);
                    % (f, 1-3\lambda)_{x_j}
                end
            end    
        case Nu % CR element
            % to-do
            return;
    end
end

[AD,BD,f,g,u,p,ufreeDof,pDof] = getbdStokesCR;


%% Record assembeling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(ufreeDof), return; end
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if length(f)+length(g) <= 5e4  % Direct solver for small size systems
        option.solver = 'direct';
    else     % Multigrid-type  solver for large size systems
        option.solver = 'diag';
    end
end
solver = option.solver;

%% Solver
switch solver
    case 'direct'
        t = cputime;
        bigA = [AD, BD'; ...
                BD, sparse(Np,Np)];
        bigF = [f; g];
        bigu = [u; p];
        bigFreeDof = [ufreeDof; pDof];
        bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
        u = bigu(1:3*Nu);
        p = bigu(3*Nu+1:end);
        residual = norm(bigF - bigA*bigu);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);        
    case 'diag'
        C = sparse(Np,Np);
%         option.freeDof = find(ufreeDof);
        option.tol = 1e-6;
        option.printlevel = 1;
        if nargin>=6
            HB = varargin{1};
        else
            HB = [];
        end        
        [u,p,info] = diagpreStokes(AD,BD,C,f,g,elem,option,HB);                                               
end

%% Post-process
if sum(pDof) ~= Np % p is unique up to a constant
    % impose the condition int(p)=0
    p = p - sum(p.*volume)/sum(volume);
end

%% Output
soln = struct('u',u,'p',p);
eqn = struct('A',AD,'B',BD,'Lap',A,'f',f,'g',g,...
             'face',face,'ufreeDof',ufreeDof,'pDof',pDof);
info.assembleTime = assembleTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdStokesCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,BD,f,g,u,p,ufreeDof,pDof] = getbdStokesCR
        %% Boundary condition of Stokes equation: CR elements
        
        %% Initial set up
        f = fxyz(:);
        g = zeros(Np,1);
        u = zeros(3*Nu,1);
        p = zeros(Np,1);
        ufreeDof = true(Nu,1); %#ok<PREALL>
        pDof = true(Np,1);
        
        if ~exist('bdFlag','var'), bdFlag = []; end
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        if ~isfield(pde,'g_R'), pde.g_R = []; end
        
        %% Part 1: Find Dirichlet dof and modify the matrix
        % Find Dirichlet boundary dof: fixedDof and pDof
        isFixedDof = false(Nu,1);
        if ~isempty(bdFlag)       % case: bdFlag is not empty
            isFixedDof(elem2face(bdFlag(:)==1)) = true; % dof on D-edges
        elseif isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
            s = accumarray(elem2face(:), 1, [Nu 1]);
            isFixedDof = (s==1);
        elseif any(isFixedDof) % pure Neumann boundary condition
            % pde.g_N could be empty which is homogenous Neumann boundary condition
            isFixedDof(1) = true;   % eliminate the kernel by enforcing u(1) = 0;
        end
        ufreeDof=~isFixedDof;
        isFixedDofAll = repmat(isFixedDof, [3,1]);
        
        % Modify the matrix
        % Build Dirichlet boundary condition into the matrix AD by enforcing
        % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
        % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.        
        Tbd = spdiags(isFixedDofAll,0,3*Nu,3*Nu);
        T = spdiags(~isFixedDofAll,0,3*Nu,3*Nu);
        AD = T*A*T + Tbd;
        BD = B*T;
        
        %% Part 2: Find boundary edges and modify the right hand side f and g
        % Find boundary edges: Neumann and Robin
        Neumann = []; Robin = []; %#ok<*NASGU>
        if ~isempty(bdFlag)
            isNeumann(elem2face((bdFlag(:)==2)|(bdFlag(:) == 3))) = true;
            isRobin(elem2face(bdFlag(:)==3)) = true;
            Neumann   = face(isNeumann,:);
            Robin     = face(isRobin,:);
        end
        if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
            % no bdFlag, only pde.g_N or pde.g_R is given in the input
            [~,Neumann] = findboundary3(elem);
            if ~isempty(pde.g_R)
                Robin = Neumann;
            end
        end
        
        % Neumann boundary condition
        if ~isempty(pde.g_N) && ~isempty(Neumann) && ~(isnumeric(pde.g_N) && all(pde.g_N == 0))
            if ~isfield(option,'gNquadorder')
                option.gNquadorder = 3;   % default order exact for linear gN
            end
            [lambdagN,weightgN] = quadpts(option.gNquadorder);
            nQuadgN = size(lambdagN,1);
            % scaled normal of face
            normalNeuman = cross(node(Neumann(:,2),:) - node(Neumann(:,1),:),...
                node(Neumann(:,3),:) - node(Neumann(:,2),:),2);
            faceArea = 0.5*sqrt(sum(normalNeuman.^2,2));
            % update RHS
            for pp = 1:nQuadgN
                ppxyz = lambdagN(pp,1)*node(Neumann(:,1),:) ... 
                      +lambdagN(pp,2)*node(Neumann(:,2),:) ...
                      +lambdagN(pp,3)*node(Neumann(:,3),:);
                gp = pde.g_N(ppxyz);
                
                for kk = 1:3 % x,y,z components
                    fxyz(isNeumann,kk) = fxyz(isNeumann,kk) + ...
                        weightgN(pp)*faceArea.*gp(:,kk);
                end
            end
        end
        f = fxyz(:);
        % The case non-empty Neumann but g_N=[] corresponds to the zero flux
        % boundary condition on Neumann edges and no modification is needed.
        
        % Dirichlet boundary conditions
        if any(isFixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))
            bdFaceCenter = (node(face(isFixedDof,1),:) ...
                          + node(face(isFixedDof,2),:) ...
                          + node(face(isFixedDof,3),:))/3;
            uD = pde.g_D(bdFaceCenter);         % bd values at centroid of faces
            u(isFixedDofAll) = uD(:); % Dirichlet bd condition is built into u
            f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
            g = g - B*u;  % the right hand side
            g = g - mean(g); % impose the compatible condition
            f(isFixedDofAll) = u(isFixedDofAll);
        end
        % The case non-empty Dirichlet but g_D=[] corresponds to the zero Dirichlet
        % boundary condition and no modification is needed.
        
        % modfiy pressure dof for pure Dirichlet
        if isempty(Neumann)
            pDof(end) = false;
        end
        
        ufreeDof = ~isFixedDofAll;              
    end % end of function getbdStokesCR
end