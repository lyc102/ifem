function [u,sigma,eqn,info] = Poisson3RT0(node,elem,bdFlag,pde,option)
%% POISSON3RT0 Poisson equation: lowest order Raviart-Thomas element in 3-D.
%
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
      center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:) ...
               +node(elem(:,4),:))/4;
      K = pde.d(center);  % take inverse sequencil.             
   end
else
    K = [];
end

%% Data structure
elemold = elem;
[elem,bdFlag] = sortelem(elem,bdFlag);  % ascend ordering
[elem2face,face] = dof3face(elem);
NT = size(elem,1); NF = size(face,1);
[Dlambda,volume,elemSign] = gradbasis3(node,elem);

%% Assemble matrix 
Nsigma = NF; Nu = NT; Ndof = Nsigma + Nu;

% M. Mass matrix for RT0 element
M = getmassmatvec3(elem2face,volume,Dlambda,'RT0',K); % ascend ordering of loc edge

% B. divergence operator
B = icdmat(double(elem2face),elemSign*[1 -1 1 -1]); % inconsistency with the induced ordering

% C. zero matrix.
C = sparse(Nu,Nu);

A = [M B';B C];

%% Assemble right hand side.% if ~isempty(bdEdge)
fu = zeros(Nu,1);
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 2;   % default order
end
if ~isempty(pde.f)
	[lambda,weight] = quadpts3(option.fquadorder);
	nQuad = size(lambda,1);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
             + lambda(p,3)*node(elem(:,3),:) ...
			 + lambda(p,4)*node(elem(:,4),:);
		fp = pde.f(pxyz);
		fu = fu - fp*weight(p);
    end
    fu = fu.*volume;
end
clear fp
F((Nsigma+1):Ndof,1) = fu;

%% Boundary condition
if ~exist('bdFlag','var'), bdFlag = []; end
[AD,F,bigu,freeDof,isPureNeumannBC] = getbd3RT0(F);
eqn = struct('M',AD(1:NF,1:NF),'B',AD(NF+1:end,1:NF),'C',AD(NF+1:end,NF+1:end),...
             'f',F(1:NF),'g',F(NF+1:end),'freeDof',freeDof,'A',AD);

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
        sigma = bigu(1:Nsigma);
        u = bigu(Nsigma+1:end); 
        info = struct('solverTime',cputime - t,'itStep',1,'error',0,'flag',0,'stopErr',0);
    case 'none'
        sigma = zeros(Nsigma,1); u = zeros(Nu,1); info =[];        
    case 'tri'
        [sigma,u,info] = tripremixPoisson(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);    
    case 'uzawapcg'
        [sigma,u,info] = uzawapcg(eqn.M,eqn.B,eqn.C,eqn.f,eqn.g,elemold);
%     case 'mg'
%         option.freeFace = freeFace;
%         option.isPureNeumannBC = isPureNeumannBC;
%         [sigma0,u,info] = mgDarcy(eqn.M,eqn.B,eqn.f,eqn.g,elemold,option);
%         sigma = bigu(1:NE);
%         sigma(freeEdge) = sigma0;
end
if isPureNeumannBC % post process for u for pure Neumann boundary condition
    ubar = sum(u.*volume)/sum(volume);
    u = u - ubar;
end

%% Output information
info.assembleTime = assembleTime; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction getbdRT0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,F,bigu,freeDof,isPureNeumannBC] = getbd3RT0(F)
    %% GETBD3RT0 Boundary conditions for Poisson equation: RT0 element in 3D.
            
    bigu = zeros(Ndof,1);

    %% Boundary conditions
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end

    %% Set up bdFlag
    if (isempty(bdFlag)) % no bdFlag information
        if ~isempty(pde.g_N)
            bdFlag = setboundary3(node,elem,'Neumann');
        elseif ~isempty(pde.g_D)
            bdFlag = setboundary3(node,elem,'Dirichlet');
        end
        % case: bdFlag = [], pde.g_D = pde.g_N =[];
        % It is equivalent to homogenous Dirichlet boundary condition
    end

    %% Find Dirichlet and Neumann dofs         
    faceSign = ones(NF,1);
    isDirichlet = false(NF,1);
    isNeumann = false(NF,1);
    if ~isempty(bdFlag)
        % Find Dirichlet and Neumann boundary edges
        isDirichlet(elem2face(bdFlag(:) == 1)) = true;
          isNeumann(elem2face(bdFlag(:) == 2)) = true;
    % Direction of boundary faces may not be the outwards normal
    % direction of the domain due to the sortelem. faceSign is
    % introduced to record this inconsistency of the ascend ordering and
    % the induced ordering.
        faceSign = ones(NF,1);
        idx = (bdFlag(:,1) ~= 0) & (elemSign == -1);% first face
        faceSign(elem2face(idx,1)) = -1;
        idx = (bdFlag(:,2) ~= 0) & (elemSign == 1); % second face
        faceSign(elem2face(idx,2)) = -1;            
        idx = (bdFlag(:,3) ~= 0) & (elemSign == -1);% third face
        faceSign(elem2face(idx,3)) = -1;
        idx = (bdFlag(:,4) ~= 0) & (elemSign == 1);% fourth face
        faceSign(elem2face(idx,4)) = -1;
    end
    Dirichlet = face(isDirichlet,:);
    Neumann = face(isNeumann,:);
    isBdDof = false(Ndof,1); 
    isBdDof(isNeumann) = true;   % for mixed method, Neumann edges are fixed
    freeDof = find(~isBdDof);
%     isFreeFace = true(NF,1);
%     isFreeFace(isNeumann) = false;
%     freeFace = find(isFreeFace);

    %% Dirichlet boundary condition (Neumann BC in mixed form)
    %   We need only modify the rhs on dof associated with Dirichlet boundary
    %   Compute the integration of g_D on the boundary Face with middle point
    %   quadrature rule. \Phi\cdot n = |F|.
    %   int_F \Phi\cdot n g_D = g_D(F_barycenter)
    %
    if ~isempty(pde.g_D) && isnumeric(pde.g_D) && (pde.g_D==0)
        pde.g_D = [];
    end        
    if ~isempty(pde.g_D) && (~isempty(Dirichlet))
        barycenter = 1/3*(node(Dirichlet(:,1),:)+node(Dirichlet(:,2),:)+ ...
                          node(Dirichlet(:,3),:));
        F(isDirichlet) = pde.g_D(barycenter).*faceSign(isDirichlet);
    end
    clear barycenter

    %% Neumann boundary condition (Dirichlet BC in mixed form)
    if ~isempty(pde.g_N) && isnumeric(pde.g_N) && (pde.g_N==0)
        pde.g_N = [];
    end    
    if ~isempty(pde.g_N) && any(isNeumann)
        % modify the rhs to include Dirichlet boundary condition
        barycenter = 1/3*(node(Neumann(:,1),:)+node(Neumann(:,2),:)+node(Neumann(:,3),:));
        ve2 = node(Neumann(:,1),:) - node(Neumann(:,3),:);
        ve3 = node(Neumann(:,2),:) - node(Neumann(:,1),:);
        ve2Crossve3 = mycross(ve2,ve3);
        faceArea = 0.5*sqrt(sum(ve2Crossve3.^2,2));
        bigu(isNeumann) = faceArea.*pde.g_N(barycenter).*faceSign(isNeumann); % 2<1,g_N>
        F = F - A*bigu;
        F(isNeumann) = bigu(isNeumann);
    end
    %% Pure Neumann boundary condition
    isPureNeumannBC = false;
    if ~any(isDirichlet) && any(isNeumann)
        freeDof = freeDof(1:end-1);  % eliminate the kernel by enforcing u(NT) = 0;
        isBdDof(end) = true;
        isPureNeumannBC = true;
        F(NF+1:end) = F(NF+1:end) - mean(F(NF+1:end)); % normalize
    end
    
    %% Modify the matrix
    %  Build Neumann boundary condition(Dirichlet BC in mixed form) into the
    %  matrix AD by enforcing  |AD(bdNode,bdNode)=I,
    %  AD(bdNode,freeNode)=0, AD(freeNode,bdNode)=0|.
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
%% TODO: Write m-lint