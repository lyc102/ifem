function [soln,eqn,info] = PoissonP3(node,elem,bdFlag,pde,option)
%% POISSONP3 Poisson equation: P3 quadratic element.
%
% u = PoissonP3(node,elem,bdFlag,pde,option) produces the cubic
%   finite element approximation of the Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
% [soln,eqn,info] = PoissonP2(node,elem,pde,bdFlag,option)
%
% The usage is the same as <a href="matlab:help Poisson">Poisson</a>. Quadratic element on a triangle is
% summarized in <a href="matlab:ifem PoissonfemrateP3">PoissonfemrateP3</a> for detail.
% 
%   Example
%     clear all
%     node = [0,0; 1,0; 1,1; 0,1];
%     elem = [2,3,1; 4,1,3];      
%     for k = 1:2
%       [node,elem] = uniformrefine(node,elem);
%     end
%     % Homogenous Dirichlet boundary condition
%     pde.f = inline('ones(size(p,1),1)','p');
%     pde.g_D = 0;
%     u = PoissonP3(node,elem,[],pde);
%     figure(1); 
%     showresult(node,elem,u);
%
% See also Poisson, Poisson3, Poisson3P2, PoissonP2 
%
% Created by Jie Zhou. Modified by Long Chen.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Preprocess
% Input arguments
if ~exist('bdFlag','var'), bdFlag = []; end
if ~exist('option','var'), option = []; end

%% Construct Data Structure
time = cputime;  % record assembling time
[elem2dof,elem2edge,edge,bdDof]  = dofP3(elem);   
% important constants
N = size(node,1);    NT = size(elem,1);     NE = size(edge,1);
Ndof = N + 2*NE + NT;

%% Compute geometric quantities and gradient of local basis
[Dlambda,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
% Since Dphi_i*Dphi_j is four degree, so four order numerical quadrature rule is used here
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'quadorder')
    option.quadorder = 5;   % default order 2(p-1)+1
end
[lambda, w] = quadpts(option.quadorder);
nQuad = size(lambda,1);
ii = zeros(55*NT,1); jj = zeros(55*NT,1); sA = zeros(55*NT,nQuad);
% generate sparse pattern
index = 0;
for i = 1:10
    for j = i:10
        ii(index+1:index+NT) = double(elem2dof(:,i)); 
        jj(index+1:index+NT) = double(elem2dof(:,j));  
        index = index + NT;
    end
end
% compute non-zeros
for p = 1:nQuad
    % Dphi at quadrature points
    Dphip(:,:,10) = 27*(lambda(p,1)*lambda(p,2)*Dlambda(:,:,3)+lambda(p,1)*lambda(p,3)*Dlambda(:,:,2)+...
                   lambda(p,3)*lambda(p,2)*Dlambda(:,:,1));   
    Dphip(:,:,1) = (27/2*lambda(p,1)*lambda(p,1)-9*lambda(p,1)+1).*Dlambda(:,:,1);           
    Dphip(:,:,2) = (27/2*lambda(p,2)*lambda(p,2)-9*lambda(p,2)+1).*Dlambda(:,:,2); 
    Dphip(:,:,3) = (27/2*lambda(p,3)*lambda(p,3)-9*lambda(p,3)+1).*Dlambda(:,:,3);
    Dphip(:,:,4) = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda(:,:,3)+...
                   lambda(p,3)*(6*lambda(p,2)-1).*Dlambda(:,:,2));  
    Dphip(:,:,5) = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda(:,:,2)+...
                   lambda(p,2)*(6*lambda(p,3)-1).*Dlambda(:,:,3)); 
    Dphip(:,:,6) = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda(:,:,1)+...
                   lambda(p,1)*(6*lambda(p,3)-1).*Dlambda(:,:,3)); 
    Dphip(:,:,7) = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda(:,:,3)+...
                   lambda(p,3)*(6*lambda(p,1)-1).*Dlambda(:,:,1)); 
    Dphip(:,:,8) = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda(:,:,2)+...
                   lambda(p,2)*(6*lambda(p,1)-1).*Dlambda(:,:,1)); 
    Dphip(:,:,9)  = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda(:,:,1)+...
                   lambda(p,1)*(6*lambda(p,2)-1).*Dlambda(:,:,2));  
    index = 0;
    for i = 1:10
        for j = i:10
            Aij = 0;
            if isempty(pde.d) || isnumeric(pde.d)
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            else
                pxy = lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:);
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.d(pxy);
            end
             if ~isempty(pde.d) && isnumeric(pde.d) % d is piecewise constant
                 Aij = pde.d.*Aij;
             end
            Aij = Aij.*area;
            sA(index+1:index+NT,p) = Aij;
            index = index + NT;
        end
    end
end
sA = sum(sA,2);
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj sA

%% Assemble right hand side by high order quadrature rule
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 6;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,w] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    phi(:,10) = 27*lambda(:,1).*lambda(:,2).*lambda(:,3);     
    phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);           
    phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2); 
    phi(:,3) = 0.5*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
    phi(:,4) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,2)-1); 
    phi(:,5) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,3)-1); 
    phi(:,6) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,3)-1);      
    phi(:,7) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,1)-1);  
    phi(:,8) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
    phi(:,9) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,2)-1);        
    bt = zeros(NT,10);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        if isfield(pde,'f') && isnumeric(pde.f)
            fp = pde.f;        % piecewise constant       
        else
            fp = pde.f(pxy);   % function handle
        end
        for j = 1:10
            bt(:,j) = bt(:,j) + w(p)*phi(p,j)*fp;
        end
    end
    bt = bt.*repmat(area,1,10);
    b = accumarray(elem2dof(:),bt(:),[Ndof 1]); 
end

%% Boundary Conditions
[AD,b,u,freeDof,isPureNeumann] = getbdP3(b);

%% Record assembeling time
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
    else            % Multigrid-type  solver for large size systems
        option.solver = 'mg';
    end
end

solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        u(freeDof) = AD(freeDof,freeDof)\b(freeDof);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        option.x0 = u;
        option.solver = 'CG';
        option.tol = 1e-9;
        [u,info] = mg(AD,b,elem,option,edge);
    case 'amg'
        option.solver = 'CG';
        option.tol = 1e-9;
        [u(freeDof),info] = amg(AD(freeDof,freeDof),b(freeDof),option);                 
end
clear phi lambda
% post-process for pure Neumann problem
if isPureNeumann
    intguh = double(sparse(NT,1));  
    [lambda,w] = quadpts(3);  %This is P3 element.
    nQuad = size(lambda,1);
    phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);           
    phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2); 
    phi(:,3) = 0.5*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
    phi(:,4) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,2)-1); 
    phi(:,5) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,3)-1); 
    phi(:,6) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,3)-1);      
    phi(:,7) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,1)-1);  
    phi(:,8) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
    phi(:,9) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,2)-1);        
    phi(:,10) = 27*lambda(:,1).*lambda(:,2).*lambda(:,3); 
    for  p = 1:nQuad
     intguh(:) = intguh(:) + w(p)*(phi(p,1)*u(elem2dof(:,1))+phi(p,2)*u(elem2dof(:,2))+phi(p,3)*u(elem2dof(:,3))...
                           + phi(p,4)*u(elem2dof(:,4))+phi(p,5)*u(elem2dof(:,5))+phi(p,6)*u(elem2dof(:,6))...
                           + phi(p,7)*u(elem2dof(:,7))+phi(p,8)*u(elem2dof(:,8))+phi(p,9)*u(elem2dof(:,9))...
                           + phi(p,10)*u(elem2dof(:,10)));
    end                 %compute the intgrable of uh.
    intguh = intguh.*area;
    uc = sum(intguh)/sum(area);
    u = u - uc;    % normalization for pure Neumann problem
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
% subfunctions getbdP3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeDof,isPureNeumann] = getbdP3(b)
    %% Boundary conditions for Poisson equation: P3 quadratic FEM.
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
    Robin = [];
    idxR = (bdFlag(:) == 3);      %index of Robin edges in bdFlag
    if any(idxR)    
        isRobin = false(NE,1);
        isRobin(elem2edge(idxR)) = true;
        Robin = edge(isRobin,:);  % Robin edges  
    end
    if ~isempty(Robin) && ~isempty(pde.g_R) && ~(isnumeric(pde.g_R) && (pde.g_R == 0))
        if ~isfield(option,'gRquadorder')
            option.gRquadorder = 8;   % we should use six order  rule.
        end
        [lambdagR,weightgR] = quadpts1(option.gRquadorder);
        nQuadgR = size(lambdagR,1);
        % cubic bases (1--3--4--2)
        bdphi = zeros(nQuadgR,4);        
        bdphi(:,1) = 0.5*(3*lambdagR(:,1)-1).*(3*lambdagR(:,1)-2).*lambdagR(:,1); 
        bdphi(:,2) = 0.5*(3*lambdagR(:,2)-1).*(3*lambdagR(:,2)-2).*lambdagR(:,2);
        bdphi(:,3) = 9/2*lambdagR(:,1).*lambdagR(:,2).*(3*lambdagR(:,1)-1);
        bdphi(:,4) = 9/2*lambdagR(:,1).*lambdagR(:,2).*(3*lambdagR(:,2)-1);        
        % length of edge
        el = sqrt(sum((node(Robin(:,1),:) - node(Robin(:,2),:)).^2,2));
        NR = size(Robin,1);
        ss = zeros(NR,4,4);
        for pp = 1:nQuadgR
            ppxy = lambdagR(pp,1)*node(Robin(:,1),:) ...
                 + lambdagR(pp,2)*node(Robin(:,2),:);
            gRp = pde.g_R(ppxy);
            for iR = 1:4
                for jR = iR:4   % only compute half of the off-diagonal part
                    ss(:,iR,jR) = ss(:,iR,jR) + ...
                    weightgR(pp)*gRp*bdphi(pp,iR).*bdphi(pp,jR);
                end
            end
        end
        ss(:) = ss(:).*repmat(el,16,1);
        Robin(:,3) = 2*find(isRobin)+N-1;  % the third one maps to corresponding dof
        Robin(:,4) = 2*find(isRobin)+N;    % the fourth one maps to corresponding dof
        index = 0;
        for iR = 1:4
            for jR = 1:4
                iiR(index+1:index+NR) = double(Robin(:,iR)); 
                jjR(index+1:index+NR) = double(Robin(:,jR)); 
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
    fixedDof = []; freeDof = [];
    isFixedDof = false(Ndof,1); 
    if ~isempty(bdFlag)     
        isDirichlet(elem2edge(bdFlag(:)==1)) = true;
        isFixedDof(edge(isDirichlet,:)) = true;
        isFixedDof(N + 2*find(isDirichlet')-1) = true;
        isFixedDof(N + 2*find(isDirichlet')) = true;
        fixedDof = find(isFixedDof);
        freeDof = find(~isFixedDof);    
    end
    if isempty(bdFlag) && ~isempty(pde.g_D) && isempty(pde.g_N) && isempty(pde.g_R)
        fixedDof = bdDof;
        isFixedDof(fixedDof) = true;
        freeDof = find(~isFixedDof);    
    end
    isPureNeumann = false;        
    if isempty(fixedDof) && isempty(Robin)  % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedDof = 1;
        freeDof = 2:Ndof;    % eliminate the kernel by enforcing u(1) = 0;
    end
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % |AD(fixedDof,fixedDof)=I, AD(fixedDof,freeDof)=0,
    % AD(freeDof,fixedDof)=0|.
    if ~isempty(fixedDof)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedDof) = 1;
        Tbd = sparse(1:Ndof,1:Ndof,bdidx,Ndof,Ndof);
        T = sparse(1:Ndof,1:Ndof,1-bdidx,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end
    
    %% Part 2: Find boundary edges and modify the load b
    % Find boundary edges: Neumann and Robin
    Neumann = [];
    if ~isempty(bdFlag)     
        idxN = (bdFlag(:) == 2);      % all Neumann edges in bdFlag        
        Neumannidx = elem2edge(idxN | idxR); % index of Neumann and Robin edges
        % since boundary integral is also needed for Robin edges
        Neumann   = edge(Neumannidx,:);
    end
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        Neumannidx = find(bdDof>N);
        Neumann = edge(Neumannidx,:);
    end
    
    % Neumann boundary condition
    if ~isempty(pde.g_N) && ~isempty(Neumann) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 6;  
        end
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        nQuadgN = size(lambdagN,1);
        % quadratic bases (1---3---4--2)
        bdphi = zeros(nQuadgN,4);        
        bdphi(:,1) = 0.5*(3*lambdagN(:,1)-1).*(3*lambdagN(:,1)-2).*lambdagN(:,1); 
        bdphi(:,2) = 0.5*(3*lambdagN(:,2)-1).*(3*lambdagN(:,2)-2).*lambdagN(:,2);
        bdphi(:,3) = 9/2*lambdagN(:,1).*lambdagN(:,2).*(3*lambdagN(:,1)-1);
        bdphi(:,4) = 9/2*lambdagN(:,1).*lambdagN(:,2).*(3*lambdagN(:,2)-1);
        % length of edge
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        ge = zeros(size(Neumann,1),4);
        for pp = 1:nQuadgN
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNp = pde.g_N(ppxy);
            ge(:,1) = ge(:,1) + weightgN(pp)*gNp*bdphi(pp,1);
            ge(:,2) = ge(:,2) + weightgN(pp)*gNp*bdphi(pp,2);
            ge(:,3) = ge(:,3) + weightgN(pp)*gNp*bdphi(pp,3);    
            ge(:,4) = ge(:,4) + weightgN(pp)*gNp*bdphi(pp,4);
        end
        ge = ge.*repmat(el,1,4);
        b(1:N) = b(1:N) + accumarray(Neumann(:), [ge(:,1); ge(:,2)],[N,1]);
        b(N+2*Neumannidx-1) = b(N+2*Neumannidx-1) + ge(:,3);
        b(N+2*Neumannidx)   = b(N+2*Neumannidx) + ge(:,4);
    end

    % Dirichlet boundary conditions
    if ~isPureNeumann && ~isempty(fixedDof) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
        % interpolation
        idx = (fixedDof > N);  % index of edge nodes
        u(fixedDof(~idx)) = pde.g_D(node(fixedDof(~idx),:)); % bd value at vertex dofs
        % for P3,  we should divide the points of edge into two parts.        
        bdEdgeIdx = fixedDof(idx) - N;
        %  First parts, the points  * is in  1---*------2
        bdEdgeMid = node(edge(isDirichlet,1),:)+(node(edge(isDirichlet,2),:) ...
                  - node(edge(isDirichlet,1),:))/3;
        u(N + bdEdgeIdx(1:2:end)) = pde.g_D(bdEdgeMid);
      %  Second parts, the points * is in  1------*---2     
        bdEdgeMid = node(edge(isDirichlet,1),:)+2*(node(edge(isDirichlet,2),:)...
                  - node(edge(isDirichlet,1),:))/3; 
        u(N + bdEdgeIdx(2:2:end)) = pde.g_D(bdEdgeMid);
        % modify the right hand side
        b = b - A*u;
    end
    if ~isPureNeumann % non-empty Dirichlet boundary condition
        b(fixedDof) = u(fixedDof);
    end
    
    % Pure Neumann boundary condition
    if isPureNeumann
        b = b - mean(b);   % compatilbe condition: sum(b) = 0
        b(1) = 0;        
    end
    end % end of getbdP3
end % end of function PoissonP3