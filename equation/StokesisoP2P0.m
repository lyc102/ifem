function [soln,eqn,info] = StokesisoP2P0(node,elem,bdFlag,pde,option)
%% STOKESISOP2P0 Stokes equation: isoP2-P0 elements.
%
%   [soln,eqn,info] = STOKESisoP2P0(node,elem,bdFlag,pde) use piceswise linear on 
%   grid h and piecewise constant on grid H = 2*h to approximate velocity 
%   u and pressure p, repectively.
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%   with
%       Dirichlet boundary condition        u = g_D  on \Gamma_D,
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.
%
% Created by Long Chen and Ming Wang, and revised at July, 2012.
% 
% See also Poisson,StokesP1P0, StokesP2P1
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('option','var'), option = []; end

%% Refine grid for P1
% [tempvar,areaC] = gradbasis(node,elem);
nodeC = node; elemC = elem; bdFlagC = bdFlag;
[tempvar,edgeC] = dofP2(elemC);
[node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
N = size(node,1);  NT = size(elem,1); Nu = N; Np = NT/4;

t = cputime;
%% Compute geometric quantities and gradient of local basis
[Dlambda,area] = gradbasis(node,elem);
areaC = sum(reshape(area,Np,4),2);

%% Assemble stiffness matrix for Laplace operator
A = sparse(Nu,Nu);
for i = 1:3
    for j = i:3
        Aij = (Dlambda(:,1,i).*Dlambda(:,1,j) + ...
               Dlambda(:,2,i).*Dlambda(:,2,j)).*area;
        if isfield(pde,'mu') && (pde.mu~=1)
            Aij = pde.mu*Aij;
        end
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [Aij; Aij],N,N);
        end
    end
end
clear Aij
A = blkdiag(A,A);

%% Assemble the matrix for divergence operator
% Dphi is constant and phi is linear and therefore 1-pt quadrature is exact
Dx = sparse(NT,N);
Dy = sparse(NT,N);
for i = 1:3
    Dx = Dx + sparse(1:NT,elem(:,i),Dlambda(:,1,i).*area,NT,N);
    Dy = Dy + sparse(1:NT,elem(:,i),Dlambda(:,2,i).*area,NT,N);
end
Dx = repmat(speye(Np,Np),1,4)*Dx;
Dy = repmat(speye(Np,Np),1,4)*Dy;
B = [-Dx -Dy];

%% Assemble right hand side by 4-points quadrature rule
f1 = zeros(Nu,1);
f2 = zeros(Nu,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f)
    % quadrature points in the barycentric coordinate
    [lambda,weight] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    ft1 = zeros(NT,3);
    ft2 = zeros(NT,3);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        % function values at quadrature points
        fp = pde.f(pxy);
        % evaluate fp outside.
        for j = 1:3
            ft1(:,j) = ft1(:,j) + fp(:,1).*lambda(p,j)*weight(p);
            ft2(:,j) = ft2(:,j) + fp(:,2).*lambda(p,j)*weight(p);
        end
    end
    ft1 = ft1.*repmat(area,1,3);
    ft2 = ft2.*repmat(area,1,3);
    f1 = accumarray(elem(:),ft1(:),[Nu 1]);
    f2 = accumarray(elem(:),ft2(:),[Nu 1]);
end

[AD,BD,f,g,u,p,ufreeDof,pDof] = getbdStokesisoP2P0;

%% Record assembeling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(ufreeDof), return; end
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if length(f)+length(g) <= 1e3  % Direct solver for small size systems
        option.solver = 'direct';
    else          % Multigrid-type  solver for large size systems
        option.solver = 'asmg';
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
        bigFreeDof = [ufreeDof; 2*Nu+pDof];
        bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
        u = bigu(1:2*Nu);
        p = bigu(2*Nu+1:end);
        residual = norm(bigF - bigA*bigu);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);        
    case 'mg'
        option.solver  = 'WCYCLE';
        [u(ufreeDof),p,info] = mgstokes(A(ufreeDof,ufreeDof),B(:,ufreeDof),f(ufreeDof),g,...
                                        u(ufreeDof),p,elemC,ufreeDof,option);         
    case 'asmg'
        [u(ufreeDof),p,info] = asmgstokes(A(ufreeDof,ufreeDof),B(:,ufreeDof),f(ufreeDof),g,...
                                          u,p,nodeC,elemC,bdFlagC,ufreeDof,option); 
end

%% Post-process
if length(pDof) ~= Np % p is unique up to a constant
    % impose the condition int(p)=0
    c = sum(p.*areaC)/sum(areaC);
    p = p - c;
end

%% Output
soln = struct('u',u,'p',p);
eqn = struct('A',AD,'B',BD,'Lap',A,'f',f,'g',g,...
             'edge',edgeC,'ufreeDof',ufreeDof,'pDof',pDof);
info.assembleTime = assembleTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdStokesisoP2P0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,BD,f,g,u,p,ufreeDof,pDof] = getbdStokesisoP2P0
        %% Boundary condition of Stokes equation: isoP2-P0 elements
        
        %% Initial set up
        g = zeros(Np,1);
        u = zeros(2*Nu,1);
        p = zeros(Np,1);
        ufreeDof = (1:Nu)';
        pDof = (1:Np)';
        if ~exist('bdFlag','var'), bdFlag = []; end
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        if ~isfield(pde,'g_R'), pde.g_R = []; end
        
        %% Part 1: Find Dirichlet dof and modify the matrix
        % Find Dirichlet boundary dof: fixedDof and pDof
        isFixedDof = false(Nu,1);
        if ~isempty(bdFlag) % case: bdFlag is not empty
            allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
            Dirichlet = allEdge((bdFlag(:) == 1),:);
            isFixedDof(Dirichlet(:)) = true;
            fixedDof = find(isFixedDof);
            ufreeDof = find(~isFixedDof);
        end
        if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N)
            fixedDof = findboundary(elem);
            isFixedDof(fixedDof) = true;
            ufreeDof = find(~isFixedDof);
        end
        if isempty(fixedDof) % pure Neumann boundary condition
            % pde.g_N could be empty which is homogenous Neumann boundary condition
            fixedDof = 1;
            ufreeDof = 2:Nu;    % eliminate the kernel by enforcing u(1) = 0;
        end
        
        % Modify the matrix
        % Build Dirichlet boundary condition into the matrix AD by enforcing
        % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
        % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.
        bdidx = zeros(2*Nu,1);
        bdidx(fixedDof) = 1;
        bdidx(Nu+fixedDof) = 1;
        Tbd = spdiags(bdidx,0,2*Nu,2*Nu);
        T = spdiags(1-bdidx,0,2*Nu,2*Nu);
        AD = T*A*T + Tbd;
        BD = B*T;
        
        %% Part 2: Find boundary edges and modify the right hand side f and g
        % Find boundary edges: Neumann and Robin
        Neumann = []; Robin = []; %#ok<*NASGU>
        if ~isempty(bdFlag)
            allEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
            Neumann = allEdge((bdFlag(:)==2)|(bdFlag(:) == 3),:);
            Robin = allEdge((bdFlag(:) == 3),:);
        end
        if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
            % no bdFlag, only pde.g_N or pde.g_R is given in the input
            [tempvar,Neumann] = findboundary(elem);
            if ~isempty(pde.g_R)
                Robin = Neumann;
            end
        end
        
        % Neumann boundary condition
        if ~isempty(pde.g_N) && ~isempty(Neumann) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
            [lambda,w] = quadpts1(3);
            nQuad = size(lambda,1);
            ve = node(Neumann(:,1),:) - node(Neumann(:,2),:); % length of edge
            edgeLength = sqrt(sum(ve.^2,2));
            % update RHS
            gex = zeros(size(Neumann,1),2);   % x-component
            gey = zeros(size(Neumann,1),2);   % y-component
            for pp = 1:nQuad
                pxy = lambda(pp,1)*node(Neumann(:,1),:)+lambda(pp,2)*node(Neumann(:,2),:);
                gp = pde.g_N(pxy);
                gex(:,1) = gex(:,1) + w(pp)*edgeLength.*gp(:,1)*lambda(pp,1);
                gex(:,2) = gex(:,2) + w(pp)*edgeLength.*gp(:,1)*lambda(pp,2);
                gey(:,1) = gey(:,1) + w(pp)*edgeLength.*gp(:,2)*lambda(pp,1);
                gey(:,2) = gey(:,2) + w(pp)*edgeLength.*gp(:,2)*lambda(pp,2);
            end
            f1(1:N) = f1(1:N) + accumarray(Neumann(:), gex(:),[N,1]);
            f2(1:N) = f2(1:N) + accumarray(Neumann(:), gey(:),[N,1]);
        end
        f = [f1; f2];
        % The case non-empty Neumann but g_N=[] corresponds to the zero flux
        % boundary condition on Neumann edges and no modification is needed.

        % Dirichlet boundary conditions
        if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
            u1 = zeros(Nu,1);
            u2 = zeros(Nu,1);
            uD = pde.g_D(node(fixedDof,:));
            u1(fixedDof) = uD(:,1);
            u2(fixedDof) = uD(:,2);          
            u = [u1;u2];
            f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition
            g = g - B*u;  % to the right hand side
            g = g - mean(g);
            f(fixedDof)    = u1(fixedDof);
            f(fixedDof+Nu) = u2(fixedDof);
        end
        % The case non-empty Dirichlet but g_D=[] corresponds to the zero Dirichlet
        % boundary condition and no modification is needed.

        % modfiy pressure dof for pure Dirichlet
        if isempty(Neumann)
            pDof = (1:Np-1)';
        end
        
        ufreeDof = [ufreeDof; Nu+ufreeDof];        
    end % end of function getbdStokesisoP2P0
end