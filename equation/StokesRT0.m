function [soln,eqn,info] = StokesRT0(node,elem,bdFlag,pde,option)
%% STOKESRT0 Stokes equations: the lowest order Raviart-Thomas element in 2D.
%
%  [soln,eqn,info] = StokesRT0(node,elem,bdFlag,pde,option)
%  uses the lowest order of Raviart-Thomas element (RT0) to approximate
%  velocity u and piecewise constant element (P0)  to approximate pressure
%  p, repectively.
%
%  We solve the following equation:
%       - grad div u + curl rot u + grad p  = f   in \Omega   
%                                  - div u  = 0   in \Omega   
%                                        u  = g   on \Gamma   
%
% ifem StokesRT0doc
%
%  See also StokesBDM1B
%
% Based on a version by Ming Wang. Revised by Lin Zhong. Discussed with Jie
% Zhou and Long Chen. Further clean up by Long Chen.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.


if ~exist('option','var'), option = []; end

%% Data structure
elemunSort = elem;
[elem,bdFlag] = sortelem(elemunSort,bdFlag);
[elem2edge,edge] = dofedge(elem);
[Clambda,area,elemSign] = curlbasis(node,elem);
N = size(node,1); NT = size(elem,1); NE = size(edge,1); 
Nu = NE; Np = NT; Ndof = Nu + Np;

%% Assemble matrices
t = cputime; % record assemble time
% Mv: Lumped mass matrix for vertex: P1 element
vecMv = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
% invMv = spdiags(1./vecMv,0,N,N);
Mv = spdiags(vecMv,0,N,N);

% Me: Mass matrix for RT0 element
Me = getmassmatvec(elem2edge,area,Clambda,'RT0');

%invMt: the inverse of Mass matrix for P0 element
invMt = spdiags(1./area,0,NT,NT);

% B: negative divergence operator
B = -icdmat(double(elem2edge),elemSign*[1 -1 1]);

% C: curl operator
C = icdmat(double(edge),[-1 1]);

% R: weak rot operator
R = spdiags(1./vecMv,0,N,N)*C'*Me;

% Vector Laplacian
A = B'*invMt*B + R'*Mv*R;

%% Assemble right hand side
locEdge = [2,3; 1 3; 1,2]; % ascend ordering
fu = zeros(Nu,1);% the right hand side of u
g = zeros(Np,1); % the right hand side of p
if ~isfield(pde,'f') || (isfield(pde,'f') && isreal(pde.f) && all(pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order is 3
end
if isfield(pde,'f') && ~isempty(pde.f)
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ... 
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for k = 1:3
            i = locEdge(k,1); j = locEdge(k,2);
            % phi_k = lambda_iClambda_j - lambda_jClambda_i;
            phi_k = lambda(p,i)*Clambda(:,:,j)-lambda(p,j)*Clambda(:,:,i);
            rhs = dot(phi_k,fp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end
    bt = bt.*repmat(area,1,3);
    fu = accumarray(elem2edge(:),bt(:),[Nu 1]);
end
clear pxy fp bt rhs phi_k psi_k

%% Boundary condition and graddiv part of A
[u,p,ufreeDof,pDof,utbd] = getbdStokesRT0;
assembleTime = cputime - t;

%% Solve the system of linear equations
% set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 1e5  % Direct solver for small size systems
        solver = 'direct';
    else             % Multigrid-type  solver for large size systems
        solver = 'mg';
    end
else
    solver = option.solver;
end
% solve the system
% get submatrices of ufreeDof
A0 = A(ufreeDof,ufreeDof);
B0 = B(:,ufreeDof);
f0 = fu(ufreeDof);
g0 = g;
if strcmp(solver,'direct') && ~isempty(ufreeDof)
    t = cputime;
    bigA = [A0, B0'; ...
            B0, sparse(Np,Np)];
    bigF = [f0; g0];
    bigu = [u; p];
    bigFreeDof = [ufreeDof; Nu+pDof];
    bigu(bigFreeDof) = bigA(1:end-1,1:end-1)\bigF(1:end-1);
    u = bigu(1:Nu);
    p = bigu(Nu+1:end);
    info.solverTime = cputime - t;
elseif strcmp(solver,'mg')
%   option.solver = 'vcycle';
    option.solver  = 'WCYCLE';
    [u(ufreeDof),p,info] = mgstokesRT0(A0,B0,f0,g0,u,p,node,elemunSort,ufreeDof,option);
end

%% Post-process
if length(pDof)~=Np % p is unique up to a constant
    % impose the condition int(p)=0
    c = sum(p.*area)/sum(area);
    p = p - c;
end
w = R*u + utbd./vecMv;

%% Output information
info.assembleTime = assembleTime;

%% Output
soln = struct('u',u,'p',p,'w',w);
eqn = struct('A',A0,'B',B0,'Me',Me,'Mv',Mv,'f',f0,'g',g0,...
             'edge',edge,'ufreeDof',ufreeDof,'pDof',pDof);
info.assembleTime = assembleTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdStokesRT0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [u,p,ufreeDof,pDof,utbd] = getbdStokesRT0
        %% Boundary condition of Stokes equation: RT0-P0 elements
        
        % Initial set up
        utbd = zeros(N,1); % the line integral of u on the boundary (u.t, tau)|_\partial \Omega
        u = zeros(Nu,1);
        p = zeros(Np,1);
        ufreeDof = (1:Nu)';
        pDof = (1:Np-1)';
        
        if ~exist('bdFlag','var'), bdFlag = []; end
        if ~isfield(pde,'g_D'), pde.g_D = []; end
        if ~isfield(pde,'g_N'), pde.g_N = []; end
        if ~isfield(pde,'g_R'), pde.g_R = []; end
        if isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
            bdFlag = [];
        end
        
        % Find Dirichlet boundary dof: fixedDof and pDof
        isFixedDof = false(Nu,1);
        if ~isempty(bdFlag)       
            isDirichlet(elem2edge(bdFlag(:)==1)) = true;
            Dirichlet = edge(isDirichlet,:);
            isFixedDof(isDirichlet) = true;
            fixedDof = find(isFixedDof);
            ufreeDof = find(~isFixedDof);
        end

        % Set up edge sign
        % edgeSign records the inconsistency of asecond orientation and
        % induced orientation for each boundary edges
        edgeSign = ones(NE,1);
        idx = (bdFlag(:,1) ~= 0 ) & (elemSign == -1) ; % the first edge is on boundary
        edgeSign(elem2edge(idx,1)) = -1;
        idx = (bdFlag(:,2) ~= 0 ) & (elemSign ==  1) ; % the second edge is on boundary
        edgeSign(elem2edge(idx,2)) = -1;
        idx = (bdFlag(:,3) ~= 0 ) & (elemSign == -1) ; % the third edge is on boundary
        edgeSign(elem2edge(idx,3)) = -1;     

        % Compute the boundary integral
        if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
            % else no bddof or g_D = 0 (no modification needed)
            % 1. Normal component of u is imposed strongly
            if (isnumeric(pde.g_D) && length(pde.g_D) == NE)
                u(fixedDof) = pde.g_D(fixedDof);
            else
                u(fixedDof) = faceinterpolate(pde.g_D,node,edge(fixedDof,:),'RT0');
            end
            % 2. Tangential component of u is imposed weakly
            [lambdagD,wgD] = quadpts1(3);
            nQuadgD = size(lambdagD,1);
            % quadrat = cputime bases 1--3--2
            bdphi = lambdagD;
            ve = node(Dirichlet(:,2),:) - node(Dirichlet(:,1),:);
            ge = zeros(size(Dirichlet,1),2);
            int_left = zeros(size(Dirichlet,1),2);
            int_right = zeros(size(Dirichlet,1),2);
            for pp = 1:nQuadgD
                ppxy = lambdagD(pp,1)*node(Dirichlet(:,1),:) ...
                     + lambdagD(pp,2)*node(Dirichlet(:,2),:);
                gDp = pde.g_D(ppxy);
                int_left = int_left + wgD(pp)*gDp*bdphi(pp,1);
                int_right = int_right + wgD(pp)*gDp*bdphi(pp,2);
            end 
            ge(:,1) = dot(int_left,ve,2).*edgeSign(fixedDof);
            ge(:,2) = dot(int_right,ve,2).*edgeSign(fixedDof);
            utbd = accumarray(Dirichlet(:), [ge(:,1); ge(:,2)],[N,1]);
        end
        
        % Modify the right hand side
        fu = fu - A*u - Me*(C*(utbd./vecMv));
         g = g - B*u;
         g = g - mean(g);
    end
end