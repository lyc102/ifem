function [u,sigma,eqn,info] = PoissonRT0(node,elem,bdFlag,pde,option)
%% POISSONRT0 Poisson equation: lowest order RT element.
%
%  [u,sigma] = PoissonRT0(node,elem,bdFlag,pde) produces an approximation of
%  the Poisson equation 
%
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N
%
%  in the mixed formulation:
%
%  Find (\sigma , u) in H_{g_N,\Gamma_N}(div,\Omega)\times L^2(\Omega) s.t. 
%
%  (d^-1\sigma,\tau) + (div \tau, u)  = <\tau*n,g_D>_{\Gamma_D} 
%  \forall \tau in H_{0,\Gamma_N}(div,\Omega) 
%  (div \sigma, v)                =  -(f,v)  
%  \forall v in L^2(\Omega) 
%
%  where 
%  H_{g,\Gamma}(div,\Omega) = {\sigma \in H(div,\Omega); \sigma*n = g 
%  on \Gamma \subset \partial\Omega }.
%
%  The unknown sigma = d*grad(u) is approximated using the lowest order
%  Raviart-Thomas element and u by piecewise constant element (with basis 1).
%
%  [u,sigma] = PoissonRT0(node,elem,bdFlag,pde,option) specifies options
%   - option.solver
%     'direct': the built in direct solver \ (mldivide)
%     'tri':    triangular preconditioner
%     'uzawapcg': PCG for the Schur complement equation
%     'none': only assemble the matrix equation but not solve
%
%   The default setting is to use the direct solver for small size problems
%   and transforming based multigrid solvers for large size problems. 
%
%  Example
%
%    squarePoissonRT0
%
% Created by Ming Wang. Reorganized by Long Chen. Change basis for u.  
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Preprocess
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end

%% Diffusion coefficient
time = cputime;  % record assembling time
if ~isfield(pde,'d'), pde.d = []; end
if isfield(pde,'d') && ~isempty(pde.d)
   if isnumeric(pde.d)
      K = pde.d;                   % d is an array
   else                            % d is a function
      center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
      K = pde.d(center);  % take inverse sequencil.             
   end
else
    K = [];
end

%% Data structure
elemold = elem;
[elem,bdFlag] = sortelem(elem,bdFlag);  % ascend ordering of elem
[elem2edge,edge] = dofedge(elem);    
NT = size(elem,1); NE = size(edge,1);
[Dlambda,area,elemSign] = gradbasis(node,elem);

%% Assemble matrix 
Nsigma = NE; Nu = NT; Ndof = Nsigma + Nu;

% M. Mass matrix for RT0 element
M = getmassmatvec(elem2edge,area,Dlambda,'RT0',K); % ascend ordering of loc edge

% B. divergence operator
B = icdmat(double(elem2edge),elemSign*[1 -1 1]); % inconsistency with the induced ordering

% C. zero matrix.
C = sparse(Nu,Nu);

A = [M B';B C];

%% Assemble right hand side.
fu = zeros(Nu,1);
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isempty(pde.f)
	[lambda,w] = quadpts(option.fquadorder);
	nQuad = size(lambda,1);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
		fu = fu - fp*w(p);  % div u = -f;
    end
    fu = fu.*area;
end
clear fp
F((Nsigma+1):Ndof,1) = fu;

%% Boundary Conditions
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,F,bigu,freeDof,isPureNeumannBC] = getbdRT0(F);
eqn = struct('M',AD(1:NE,1:NE),'B',AD(NE+1:end,1:NE),'C',AD(NE+1:end,NE+1:end),...
             'f',F(1:NE),'g',F(NE+1:end),'freeDof',freeDof,'A',AD);

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the linear system.
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else         % MGCG  solver for large size systems
        option.solver = 'tri';
    end
% elseif strcmp(option.solver,'mg')
%     option.solver = 'tripremixPoisson';    
end
solver = option.solver;
% solve
switch lower(solver)
    case 'direct'
        t = cputime;
        bigu(freeDof) = AD(freeDof,freeDof)\F(freeDof);
        sigma = bigu(1:NE);
        u = bigu(NE+1:end); 
        info = struct('solverTime',cputime - t,'itStep',1,'error',0,'flag',0,'stopErr',0);
    case 'none'
        sigma = zeros(NE,1); u = zeros(NT,1); info =[];        
    case 'tri'
        [sigma,u,info] = tripremixPoisson(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);    
    case 'uzawapcg'
        [sigma,u,info] = uzawapcg(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);
%     case 'mg'
%         option.freeEdge = freeEdge;
%         option.isPureNeumannBC = isPureNeumannBC;
%         [sigma0,u,info] = mgDarcy(eqn.M,eqn.B,eqn.f,eqn.g,elemold,option);
%         sigma = bigu(1:NE);
%         sigma(freeEdge) = sigma0;
end
if isPureNeumannBC % post process for u for pure Neumann boundary condition
    ubar = sum(u.*area)/sum(area);
    u = u - ubar;
end

%% Output information
info.assembleTime = assembleTime; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction getbdRT0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,F,bigu,freeDof,isPureNeumannBC] = getbdRT0(F)
    %% GETBDRT0 Boundary conditions for Poisson equation: RT0 element.
    %
    %  Created by Ming Wang. Improved the check of edgeSign by Long Chen.

    bigu = zeros(Ndof,1);
    
    %% No boundary conditions
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end

    %% Set up bdFlag
    if isempty(bdFlag) % no bdFlag information
       if ~isempty(pde.g_N) % case: Neumann
           bdFlag = setboundary(node,elem,'Neumann');
       elseif ~isempty(pde.g_D) % case: Dirichlet
           bdFlag = setboundary(node,elem,'Dirichlet');
       end
    end

    %% Find Dirichlet and Neumann dofs 
    edgeSign = ones(NE,1);
    isDirichlet = false(NE,1);
    isNeumann = false(NE,1);
    if ~isempty(bdFlag)
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
          isNeumann(elem2edge(bdFlag(:)==2)) = true;
        % Direction of boundary edges may not be the outwards normal direction 
        % of the domain. edgeSign is introduced to record this inconsistency.
        edgeSign = ones(NE,1);
        idx = (bdFlag(:,1) ~= 0) & (elemSign == -1);% first edge is on boundary
        edgeSign(elem2edge(idx,1)) = -1;
        idx = (bdFlag(:,2) ~= 0) & (elemSign == 1); % second edge is on boundary
        edgeSign(elem2edge(idx,2)) = -1;            
        idx = (bdFlag(:,3) ~= 0) & (elemSign == -1);% third edge is on boundary
        edgeSign(elem2edge(idx,3)) = -1;
    end
    Dirichlet = edge(isDirichlet,:);
    Neumann = edge(isNeumann,:); 
    isBdDof = false(Ndof,1); 
    isBdDof(isNeumann) = true;   % for mixed method, Neumann edges are fixed
    freeDof = find(~isBdDof);
%     isFreeEdge = true(NE,1);
%     isFreeEdge(isNeumann) = false;
%     freeEdge = find(isFreeEdge);
    
    %% Dirichlet boundary condition (Neumann BC in the mixed form)
    %   We need only modify the rhs on dof associated with Dirichlet
    %   boundary. Compute the int_e g_D \phi_e\cdot n ds on the boundary.
    if ~isempty(pde.g_D) && isnumeric(pde.g_D) && (pde.g_D==0)
        pde.g_D = [];
    end
    if ~isempty(pde.g_D) && any(isDirichlet) 
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for quadratic gN
        end
        [lambda,w] = quadpts1(option.gNquadorder);
        nQuad = size(lambda,1);
        % <\phi\cdot n, g_D> = 1/|e_{i,j}|\int e_{i,j} g_D ds        
        for ip = 1:nQuad
        	pxy = lambda(ip,1)*node(Dirichlet(:,1),:)+...
                  lambda(ip,2)*node(Dirichlet(:,2),:);               
            F(isDirichlet) = F(isDirichlet) + w(ip)*pde.g_D(pxy);
        end
        F(isDirichlet) = F(isDirichlet).*edgeSign(isDirichlet);
        % no edge length since the basis of sigma contains it.
    end

    %% Neumann boundary condition (Dirichlet BC in mixed form)
    % We compute the integral int_e g_N ds and assign to boundary dof
    if ~isempty(pde.g_N) && isnumeric(pde.g_N) && (pde.g_N==0)
        pde.g_N = [];
    end    
    if ~isempty(pde.g_N) && any(isNeumann)
        % modify the rhs to include Dirichlet boundary condition 
        ve = node(Neumann(:,1),:)-node(Neumann(:,2),:);
        edgeLength = sqrt(sum(ve.^2,2)); 
        % compute the integral int_e g_N ds
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for quadratic gN
        end
        [lambda,w] = quadpts1(option.gNquadorder);
        nQuad = size(lambda,1);
        bigu(isNeumann) = 0;
        for ip = 1:nQuad
        	pxy = lambda(ip,1)*node(Neumann(:,1),:)+...
                  lambda(ip,2)*node(Neumann(:,2),:);               
            bigu(isNeumann) = bigu(isNeumann) + w(ip)*pde.g_N(pxy).*edgeLength;
        end
        % correct the sign
        bigu(isNeumann) = bigu(isNeumann).*edgeSign(isNeumann);
        F = F - A*bigu;
        F(isNeumann) = bigu(isNeumann);
    end
    
    %% Pure Neumann boundary condition
    isPureNeumannBC = false;
    if ~any(isDirichlet) && any(isNeumann)
        freeDof = freeDof(1:end-1);  % eliminate the kernel by enforcing u(NT) = 0;
        isBdDof(end) = true;
        isPureNeumannBC = true;
        F(NE+1:end) = F(NE+1:end) - mean(F(NE+1:end)); % normalize
%         F(end) = 0;
    end

    %% Modify the matrix
    %  Build Neumann boundary condition(Dirichlet BC in mixed form) into the
    %  matrix AD by enforcing  |AD(bdNode,bdNode)=I, 
    %  AD(bdNode,FreeNode)=0, AD(FreeNode,bdNode)=0|.
    if any(isBdDof)
       bdidx = zeros(Ndof,1); 
       bdidx(isBdDof) = 1;
       Tbd = spdiags(bdidx,0,Ndof,Ndof);
       T = spdiags(1-bdidx,0,Ndof,Ndof);
       AD = T*A*T + Tbd;
    else
       AD = A;
    end
    end % end of getbdRT0
end