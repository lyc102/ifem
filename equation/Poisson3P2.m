function [soln,eqn,info] = Poisson3P2(node,elem,bdFlag,pde,option,varargin)
%% POISSON3P2 Poisson equation: P2 quadratic element in 3-D
%
% u = Poisson3P2(node,elem,bdFlag,pde,option) produces the quadratic
%   finite element approximation of the Poisson equation in 3-D
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
% [soln,eqn,info] = Poisson3P2(node,elem,pde,bdFlag,option)
%
% The usage is the same as <a href="matlab:help Poisson">Poisson</a>. Quadratic element on a tetrahedron is
% summarized in <a href="matlab:ifem Poisson3P2femrate">Poisson3P2femrate</a> for detail.
%
%   Example
%
%     cubePoissonP2;
%
%     Poisson3P2femrate;
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Parameters
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end
N = size(node,1); 
NT = size(elem,1); 

%% Construct Data Structure
time = cputime;  % record assembling time
[elem2dof,edge] = dof3P2(elem);
Ndof = N + size(edge,1);

%% Compute geometric quantities and gradient of local basis
[Dlambda,volume] = gradbasis3(node,elem);

%% Assemble stiffness matrix
% generate sparse pattern
ii = zeros(55*NT,1); 
jj = zeros(55*NT,1); 
index = 0;
for i = 1:10
    for j = i:10
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));  
        index = index + NT;
    end
end
% quadrature points
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'quadorder')
    % diffusion is piecewise constant
    option.quadorder = 2;        % default order
    if ~isempty(pde.d) && isnumeric(pde.d) % numerical diffusion
        option.quadorder = 3;    % exact for linear diffusion coefficient
    end
end
[lambda, w] = quadpts3(option.quadorder);
nQuad = size(lambda,1);
% compute non-zeros
sA = zeros(55*NT,nQuad);
for p = 1:nQuad
    % Dphi at quadrature points
    Dphip(:,:,10) = 4*(lambda(p,3)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,3));
    Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
    Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
    Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
    Dphip(:,:,4) = (4*lambda(p,4)-1).*Dlambda(:,:,4);
    Dphip(:,:,5) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,1));
    Dphip(:,:,7) = 4*(lambda(p,1)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,1));
    Dphip(:,:,8) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
    Dphip(:,:,9) = 4*(lambda(p,2)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,2));
    index = 0;
    for i = 1:10
        for j = i:10
            Aij = 0;
            if isempty(pde.d) || isnumeric(pde.d)
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            else
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ...
                     + lambda(p,4)*node(elem(:,4),:);
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.d(pxyz);                
            end
            if ~isempty(pde.d) && isnumeric(pde.d) % d is piecewise constant
                Aij = pde.d.*Aij;
            end
            Aij = Aij.*volume;
            sA(index+1:index+NT,p) = Aij;
            index = index + NT;
        end
    end
end
sA = sum(sA,2);
% assemble the matrix
diagIdx = (ii == jj);   
upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj sA

%% Assemble right hand side by quadrature rule
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 4;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,w] = quadpts3(option.fquadorder);
    nQuad = size(lambda,1);
    phi(:,10)= 4*lambda(:,3).*lambda(:,4);
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phi(:,4) = lambda(:,4).*(2*lambda(:,4)-1);
    phi(:,5) = 4*lambda(:,1).*lambda(:,2);
    phi(:,6) = 4*lambda(:,1).*lambda(:,3);
    phi(:,7) = 4*lambda(:,1).*lambda(:,4);
    phi(:,8) = 4*lambda(:,2).*lambda(:,3);
    phi(:,9) = 4*lambda(:,2).*lambda(:,4);
    bt = zeros(NT,10);  
    for p = 1:nQuad
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:) ...
             + lambda(p,4)*node(elem(:,4),:);
        if isfield(pde,'f') && isnumeric(pde.f)
            fp = pde.f;        % piecewise constant       
        else
            fp = pde.f(pxyz);   % function handle
        end
        for j = 1:10
            bt(:,j) = bt(:,j) + w(p)*phi(p,j)*fp;
        end
    end
    bt = bt.*repmat(volume,1,10);
    b = accumarray(elem2dof(:),bt(:),[Ndof 1]);
end


%% Boundary Conditions
[AD,b,u,freeDof,isPureNeumann] = getbd3P2(b);

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeDof), return; end
% Set up solver
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
if isPureNeumann
    option.solver = 'mg';
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        t = cputime;
        u(freeDof) = AD(freeDof,freeDof)\b(freeDof);
        residual = norm(b - AD*u);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        option.x0 = u;
        option.solver = 'CG';
        if nargin>=6
            HB = varargin{1};
        else
            HB = [];
        end
        [u,info] = mg(AD,b,elem,option,HB); 
    case 'amg'
        option.solver = 'CG';
        [u(freeDof),info] = amg(AD(freeDof,freeDof),b(freeDof),option);                 
end
% post-process for pure Neumann problem
% still questionable for quadratic P2 in 3D
if isPureNeumann
    intu = (sum(u(elem2dof(:,5:10)),2)/5 - sum(u(elem2dof(:,1:4)),2)/20).*volume;
    uc = sum(intu)/sum(volume);
    u = u - uc;   % int u = 0
end

%% Compute Du
Du = [];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'edge',edge,'freeDof',freeDof,'Lap',A);
    info.assembleTime = assembleTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdP2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeDof,isPureNeumann] = getbd3P2(b)
%% Boundary conditions for Poisson equation: P2 quadratic FEM in 3-D.
    %
    % The set up of boundary condition consists of two parts: 
    %
    % 1) Modify the matrix for Dirichlet boundary nodes, which are not degree
    % of freedom. Values at these nodes are evaluatation of pde.g_D. The
    % original stiffness matrix A is turn into the matrix AD by enforcing
    % AD(fixedDof,fixedDof)=I, AD(fixedDof,freeDof)=0, AD(freeDof,fixedDof)=0.
    %
    % 2) Modify the right hand side b. The Neumann boundary integral is added
    % to b. For Dirichlet boundary ndoes, b(fixedDof) is the evaluation of
    % pde.g_D.
    %
    % Special attentation should be given for the pure Neumann boundary
    % condition. To enforce the compatible condition, the vector b should have
    % mean value zero. To avoid a singular matrix, the 1st node is chosen as
    % fixedDof. 
    %
    % The order of assigning Neumann and Dirichlet boundary condition is
    % important to get the right setting at the intersection nodes of Dirichlet
    % and Neumann boundary edges.

    u = zeros(Ndof,1);    
    %% Set up boundary and basic parameter
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end
    
    %% Part 1: Modify the matrix for Dirichlet and Robin condition
    % Robin boundary condition
    % Find Robin face dofs
    RobinFace2dof = [];
    if ~isempty(bdFlag) 
        RobinFace2dof = [elem2dof(bdFlag(:,1) == 3,[2,3,4,10,9,8]); ...  
                         elem2dof(bdFlag(:,2) == 3,[1,4,3,10,6,7]); ...
                         elem2dof(bdFlag(:,3) == 3,[1,2,4,9,7,5]); ...
                         elem2dof(bdFlag(:,4) == 3,[1,3,2,8,5,6])];
    end
    % Assemble the mass matrix for Robin boundary condition
    if ~isempty(RobinFace2dof) && ~isempty(pde.g_R)
        v12 = node(RobinFace2dof(:,2),:)-node(RobinFace2dof(:,1),:);
        v13 = node(RobinFace2dof(:,3),:)-node(RobinFace2dof(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        if ~isfield(option,'gRquadorder')
            option.gRquadorder = 4; 
        end
        [lambdagR,weightgR] = quadpts(option.gRquadorder);
        bdphi(:,6) = 4*lambdagR(:,1).*lambdagR(:,2);
        bdphi(:,1) = lambdagR(:,1).*(2*lambdagR(:,1)-1);
        bdphi(:,2) = lambdagR(:,2).*(2*lambdagR(:,2)-1);
        bdphi(:,3) = lambdagR(:,3).*(2*lambdagR(:,3)-1);
        bdphi(:,4) = 4*lambdagR(:,2).*lambdagR(:,3);
        bdphi(:,5) = 4*lambdagR(:,3).*lambdagR(:,1);
        nQuadgR = size(lambdagR,1);
        NR = size(RobinFace2dof,1);
        ss = zeros(NR,6,6);
        % int g_R phi_i phi_j
        for pp = 1:nQuadgR
            % quadrature points in the x-y coordinate
            ppxyz = lambdagR(pp,1)*node(RobinFace2dof(:,1),:) ...
                  + lambdagR(pp,2)*node(RobinFace2dof(:,2),:) ...
                  + lambdagR(pp,3)*node(RobinFace2dof(:,3),:);
            gRp = pde.g_R(ppxyz);
            for iR = 1:6
                for jR = iR:6
                    ss(:,iR,jR) = ss(:,iR,jR) + ...
                    weightgR(pp)*gRp*bdphi(pp,iR).*bdphi(pp,jR);
                end
            end
        end
        ss(:) = ss(:).*repmat(area,36,1);     
        % assemble
        index = 0;
        for iR = 1:6
            for jR = 1:6
                iiR(index+1:index+NR) = double(RobinFace2dof(:,iR)); 
                jjR(index+1:index+NR) = double(RobinFace2dof(:,jR)); 
                if jR>=iR
                    ssR(index+1:index+NR) = ss(:,iR,jR);
                else
                    ssR(index+1:index+NR) = ss(:,jR,iR);
                end
                index = index + NR;
            end
        end
        A = A + sparse(iiR,jjR,ssR,Ndof,Ndof);        
    end
               
    % Find Dirichlet boundary dof: fixedDof
    fixedDof = []; 
    freeDof = [];
    isFixedDof = false(Ndof,1); 
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        % no bdFlag, only pde.g_D is given. Then set as Dirichlet boundary condition only
        bdFlag = setboundary3(node,elem,'Dirichlet');        
    end
    if ~isempty(bdFlag)
        % Find boundary edges and nodes
        isFixedDof(elem2dof(bdFlag(:,1) == 1,[2,3,4,8,9,10])) = true;
        isFixedDof(elem2dof(bdFlag(:,2) == 1,[1,3,4,6,7,10])) = true;
        isFixedDof(elem2dof(bdFlag(:,3) == 1,[1,2,4,5,7,9])) = true;
        isFixedDof(elem2dof(bdFlag(:,4) == 1,[1,2,3,5,6,8])) = true;
        fixedDof = find(isFixedDof);
        freeDof = ~isFixedDof;    
    end
    
    % Modify the matrix for different boundary conditions
    AD = A;
    % Dirichlet boundary condition
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0.
    if ~isempty(fixedDof)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedDof) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    end
    % Neumann boundary condition    
    isPureNeumann = false;        
    if isempty(fixedDof) && isempty(RobinFace2dof)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        AD(1,1) = AD(1,1) + 1e-6; % eliminate the kernel        
%         fixedDof = 1;
%         freeDof = 2:Ndof;    % eliminate the kernel by enforcing u(1) = 0;
    end

    %% Part 2: Find boundary faces and modify the load b
    % Find boundary faces bdFace for Neumann boundary condition
    bdFace2dof = [];
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        bdFlag = setboundary3(node,elem,'Neumann');        
    end
    if ~isempty(bdFlag)
        % Find boundary faces for both Neumann and Robin conditons
        bdFace2dof = [elem2dof(bdFlag(:,1) >= 2,[2,3,4,10,9,8]); ...  
                      elem2dof(bdFlag(:,2) >= 2,[1,4,3,10,6,7]); ...
                      elem2dof(bdFlag(:,3) >= 2,[1,2,4,9,7,5]); ...
                      elem2dof(bdFlag(:,4) >= 2,[1,3,2,8,5,6])];
    end
    % Neumann boundary condition
    if ~isempty(bdFace2dof) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        v12 = node(bdFace2dof(:,2),:)-node(bdFace2dof(:,1),:);
        v13 = node(bdFace2dof(:,3),:)-node(bdFace2dof(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 3; 
        end
        [lambdagN,weightgN] = quadpts(option.gNquadorder);
        phigN(:,6) = 4*lambdagN(:,1).*lambdagN(:,2);
        phigN(:,1) = lambdagN(:,1).*(2*lambdagN(:,1)-1);
        phigN(:,2) = lambdagN(:,2).*(2*lambdagN(:,2)-1);
        phigN(:,3) = lambdagN(:,3).*(2*lambdagN(:,3)-1);
        phigN(:,4) = 4*lambdagN(:,2).*lambdagN(:,3);
        phigN(:,5) = 4*lambdagN(:,3).*lambdagN(:,1);
        nQuadgN = size(lambdagN,1);
        gf = zeros(size(bdFace2dof,1),6);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxyz = lambdagN(pp,1)*node(bdFace2dof(:,1),:) ...
                  + lambdagN(pp,2)*node(bdFace2dof(:,2),:) ...
                  + lambdagN(pp,3)*node(bdFace2dof(:,3),:);
            gNp = pde.g_N(ppxyz);
            for iN = 1:6
                gf(:,iN) = gf(:,iN) + weightgN(pp)*phigN(pp,iN)*gNp;
            end
        end
        gf = gf.*repmat(area,1,6);
        b = b + accumarray(bdFace2dof(:),gf(:),[Ndof,1]); 
    end
    % The case with non-empty Neumann edges but g_N=0 or g_N=[] corresponds to
    % the zero flux boundary condition on Neumann edges and no modification of
    % A,u,b is needed.

    % Dirichlet boundary condition
    if ~isPureNeumann && ~isempty(fixedDof) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
        idx = (fixedDof > N);                         % index of edge nodes
        u(fixedDof(~idx)) = pde.g_D(node(fixedDof(~idx),:)); % bd value at vertex dofs
        bdEdgeIdx = fixedDof(idx) - N;
        bdEdgeMid = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2;
        u(fixedDof(idx)) = pde.g_D(bdEdgeMid);
        b = b - A*u;
    end
    if ~isPureNeumann && ~isempty(fixedDof) % non-empty Dirichlet boundary condition
        b(fixedDof) = u(fixedDof);
    end
    % The case with non-empty Dirichlet nodes but g_D=0 or g_D=[] corresponds
    % to the zero Dirichlet boundary condition and no modification of u,b is
    % needed.

    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
%         b(1) = 0;          % 1 is fixedDof and set u(1) = 0
    end
    end % end of getbd3P2
end % end of function Poisson3P2