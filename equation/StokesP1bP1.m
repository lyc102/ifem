function [soln,eqn,info] = StokesP1bP1(node,elem,bdFlag,pde,option)
%% STOKESMINI Stokes equation: P1bP1 elements.
%
%   [soln,eqn,info] = STOKESP1bP1(node,elem,bdFlag,pde) use continous P1+bubble
%   function and continous P1 element to approximate velocity u and 
%   pressure p, repectively.
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0 in \Omega,
%   with
%       Dirichlet boundary condition        u = g_D  on \Gamma_D,
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.
%
%  Direct discretization of P1bP1 element yields the following system:
%
%             | A  B' ||u |  =  |f|
%             | B  0  ||p |  =  |0|
%
%  The only difference between StokesP1bP1 and StokesMini is that the whole
%  system is assembled, instead of eliminating the bubble DOF in StokesMini. 
%  This system is better used for Multigrid. 
%
%  Created by Ming Wang and Long Chen at Jan, 2012.
%

if ~exist('option','var'), option = []; end

%% Construct Data Structure
[tempvar,edge] = dofedge(elem);
N = size(node,1);  NT = size(elem,1); Nu = N+NT; Np = N;

t = cputime;
%% Compute geometric quantities and gradient of locA basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix for Laplace
%** Linear element part
A = sparse(N,N);
Dij = zeros(NT,1);
for i = 1:3
    for j = i:3
        Aij = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if (j==i)
            Dij = Dij + Aij;
            A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],...
                [Aij; Aij],N,N);
        end
    end
end
%** Buble function part: 
%      D = sum_i(\nabla \lambda_i,\nabla \lambda_i)/180
D = spdiags(Dij/180,0,NT,NT)*27^2;
A = blkdiag(A,D,A,D);

clear Aij Dij

%% Assemble the matrix for divergence operator
Dx = sparse(N,N); Dy = sparse(N,N); Dxb = sparse(N,NT); Dyb = sparse(N,NT);
for i = 1:3 % loop for p index
    % Linear part
    for j = 1:3 % loop for u index
        Dx = Dx + sparse(elem(:,i), elem(:,j), 1/3*Dphi(:,1,j).*area, N, N);
        Dy = Dy + sparse(elem(:,i), elem(:,j), 1/3*Dphi(:,2,j).*area, N, N);
    end
    % Buble part:  \int div (phi_b, 0)\lambda_j =
    %             -\int (\lambda_j)_x phi_b = -1/60*\int (\lambda_j)_x
    Dxb = Dxb + sparse(elem(:,i), (1:NT), -1/60*Dphi(:,1,i).*area, N, NT)*27;
    Dyb = Dyb + sparse(elem(:,i), (1:NT), -1/60*Dphi(:,2,i).*area, N, NT)*27;
end
B = -[Dx Dxb Dy Dyb];
clear Dx Dy Dxb Dyb Dphi

%% Assemble right hand side
f1 = zeros(N,1);
f2 = zeros(N,1);
g1 = zeros(NT,1);
g2 = zeros(NT,1);
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f)
    % Linear part: by 4-points quadrature rule
    mid1 = (node(elem(:,2),:) + node(elem(:,3),:))/2;
    mid2 = (node(elem(:,3),:) + node(elem(:,1),:))/2;
    mid3 = (node(elem(:,1),:) + node(elem(:,2),:))/2;
    fmid1 = pde.f(mid1);    fmid2 = pde.f(mid2);    fmid3 = pde.f(mid3);
    ft1(:,1) = area.*(fmid2(:,1) + fmid3(:,1))/6;
    ft1(:,2) = area.*(fmid3(:,1) + fmid1(:,1))/6;
    ft1(:,3) = area.*(fmid1(:,1) + fmid2(:,1))/6;
    ft2(:,1) = area.*(fmid2(:,2) + fmid3(:,2))/6;
    ft2(:,2) = area.*(fmid3(:,2) + fmid1(:,2))/6;
    ft2(:,3) = area.*(fmid1(:,2) + fmid2(:,2))/6;
    f1 = accumarray(elem(:),ft1(:),[N 1]);
    f2 = accumarray(elem(:),ft2(:),[N 1]);
    % Buble part
    [lambda,weight] = quadpts(4);
    phi = lambda(:,1).*lambda(:,2).*lambda(:,3)*27;
    nQuad = size(lambda,1);
    g1 = zeros(NT,1); g2 = zeros(NT,1);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        % function values at quadrature points
        fp = pde.f(pxy);
        g1 = g1+fp(:,1).*phi(p)*weight(p).*area;
        g2 = g2+fp(:,2).*phi(p)*weight(p).*area;
    end
end

%% Boundary condition
[AD,BD,f,g,u,p,ufreeDof,pDof] = getbdStokesP1bP1;

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
                                        u(ufreeDof),p,elem,ufreeDof,option);         
    case 'asmg'
        [u(ufreeDof),p,info] = asmgstokes(A(ufreeDof,ufreeDof),B(:,ufreeDof),f(ufreeDof),g,...
                                          u,p,node,elem,bdFlag,ufreeDof,option); 
end

%% Post-process
if length(pDof)~=Np % p is unique up to a constant
    % impose the condition int(p)=0
    c = sum(mean(p(elem),2).*area)/sum(area);
    p = p - c;
end

%% Output
soln = struct('u',u,'p',p);
eqn = struct('A',AD,'B',BD,'Lap',A,'f',f,'g',g,...
             'edge',edge,'ufreeDof',ufreeDof,'pDof',pDof);
info.assembleTime = assembleTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdStokesP1bP1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,BD,f,g,u,p,ufreeDof,pDof] = getbdStokesP1bP1
        %% Boundary condition of Stokes equation: P1bP1 element
        
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
        bdidx([fixedDof; Nu+fixedDof]) = 1;
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
        
        f = [f1; g1; f2; g2];
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
    end % end of function getbdStokesP1bP1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end