function [sigma,u,eqn,info] = HodgeLaplacian3E(node,elem,bdFlag,pde,option)
%% HODEGELAPLACIAN3E Hodge Laplacian of edge element in three dimensions
%
% [sigma,u] = HodgeLaplacian3E(node,elem,bdFlag,pde) computes the mixed
% finite element method approximation of the Hodge Laplacian equation
%
% Multigrid solvers are added by Jie Zhou.


if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

tic;
%% Data structure
elemMG  = elem;     % save elem and bdFlag for multigrid
bdFlagMG = bdFlag;
[elem,bdFlag] = sortelem3(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dof3edge(elem);
[Dlambda,volume,elemSign] = gradbasis3(node,elem);
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
Nsigma = N; Nu = NE; Ndof = Nsigma + Nu;

%% Assemble matrix 
% Mass matrices
Mv = getmassmat3(node,elem,volume);
Me = getmassmatvec3(elem2edge,volume,Dlambda,'ND0');
% G. gradient operator: P1 -> (ND1)'
G = Me*icdmat(double(edge),[-1 1]);  % gradient matrix
% R'MR: curlcurl operator
curlPhi(:,:,6) = 2*mycross(Dlambda(:,:,3),Dlambda(:,:,4),2);
curlPhi(:,:,1) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,2),2);
curlPhi(:,:,2) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,3),2);
curlPhi(:,:,3) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,4),2);
curlPhi(:,:,4) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,3),2);
curlPhi(:,:,5) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,4),2);
ii = zeros(21*NT,1); jj = zeros(21*NT,1); sC = zeros(21*NT,1); 
index = 0;
for i = 1:6
    for j = i:6
        % local to global index map
        % curl-curl matrix
        Cij = dot(curlPhi(:,:,i),curlPhi(:,:,j),2).*volume;
        ii(index+1:index+NT) = double(elem2edge(:,i)); 
        jj(index+1:index+NT) = double(elem2edge(:,j));
        sC(index+1:index+NT) = Cij;
        index = index + NT;
    end
end
clear curlPhi % clear large size array
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
C = sparse(ii(diagIdx),jj(diagIdx),sC(diagIdx),NE,NE);
CU = sparse(ii(upperIdx),jj(upperIdx),sC(upperIdx),NE,NE);
C = C + CU + CU';
% big matrix for the mixed system
A = [-Mv G'; G C];

%% Assemble right hand side
F = zeros(Ndof,1);
if ~isfield(pde,'f') || (isfield(pde,'f') && isreal(pde.f) && all(pde.f==0))
    pde.f = [];
end
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order is 3
end
if ~isempty(pde.f)
    [lambda,w] = quadpts3(option.fquadorder);
    nQuad = size(lambda,1);
    bt = zeros(NT,6);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ... 
             + lambda(p,3)*node(elem(:,3),:) ... 
             + lambda(p,4)*node(elem(:,4),:);
        fp = pde.f(pxyz);
%         locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        for k = 1:6
            i = locEdge(k,1); j = locEdge(k,2);
            % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
            phi_k = lambda(p,i)*Dlambda(:,:,j)-lambda(p,j)*Dlambda(:,:,i);
            rhs = dot(phi_k,fp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end
    bt = bt.*repmat(volume,1,6);
    F = accumarray(N+elem2edge(:),bt(:),[Ndof 1]);
end
clear pxyz fp bt rhs phi_k Dlambda

%% Boundary Conditions
[AD,F,u,sigma,freeNode,freeEdge,freeDof,isPureNeumann] = getbdHodgeLap3E(F);
% reduce to free dofs
G  = G(freeEdge,freeNode);
C  = C(freeEdge,freeEdge);
Mv = Mv(freeNode,freeNode);
f = F(freeNode);
g = F(N+freeEdge);

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

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
        tic;
        temp = zeros(Ndof,1);
        temp(freeDof) = A(freeDof,freeDof)\F(freeDof);
        sigma(freeNode) = temp(freeNode);
        u(freeEdge) = temp(freeEdge+N);
        residual = norm(F - AD*[sigma; u]);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'TRI' % GMRES with a triangular preconditioner 
        [sigma0,u0,info] =  tripreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
    case 'SDIAG' % GMRES with a signed diagonal preconditioner 
        [sigma0,u0,info] =  sdiagpreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
    case 'APPF' % GMRES with a signed diagonal preconditioner 
        [sigma0,u0,info] =  appfactpreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
    case 'DIAG' % MINRES with a diagonal preconditioner
        [sigma0,u0,info] =  diapreHodgeLapE(Mv,G,C,f,g,node,elemMG,bdFlagMG,option.mg);
        sigma(freeNode)  = sigma0;
        u(freeEdge)      = u0;
end
if isPureNeumann
   u = u - mean(u); % normalize u 
end

%% Output
eqn = struct('A',AD,'Mv',Mv,'Me',Me,'G',G,'C',C,'f',F,'edge',edge,...
             'freeNode',freeNode,'freeEdge',freeEdge);
info.assembleTime = assembleTime;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdHodgeLap3E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,f,u,sigma,freeNode,freeEdge,freeDof,isPureNeumann] = getbdHodgeLap3E(f)

    u = zeros(NE,1);  sigma = zeros(N,1);
    %% Neumann boundary condition
    [fixedNode,Neumann,isBdNode] = findboundary3(elem,bdFlag);     %#ok<*ASGLU>
    % Part u: int_F (curl u cross n, v)
    g = zeros(Ndof,1);
    if ~isfield(pde,'gcurlu') || (isnumeric(pde.gcurlu) && all(pde.gcurlu == 0))
        pde.gcurlu = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.gcurlu)
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
    % Part sigma: int_F (u dot n, tau)
    if ~isfield(pde,'gun') || (isnumeric(pde.gun) && all(pde.gun == 0))
        pde.gun = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.gun)
        v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
        v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        % three middle points rule
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 3;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts(option.gNquadorder);
        phigN = lambdagN;                 % linear bases
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),3);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxyz = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                  + lambdagN(pp,2)*node(Neumann(:,2),:) ...
                  + lambdagN(pp,3)*node(Neumann(:,3),:);
            gNp = pde.gun(ppxyz);
            for iN = 1:3
                ge(:,iN) = ge(:,iN) + weightgN(pp)*phigN(pp,iN)*gNp;
            end
        end
        ge = ge.*repmat(area,1,3);
        f = f + accumarray(Neumann(:), ge(:),[Ndof,1]); 
    end

    %% Dirichlet boundary condition
    % Find Dirichlet boundary dof: fixedDof
    if ~isempty(bdFlag)
        % Find boundary edges and nodes
        isFixedEdge = false(NE,1);
        isFixedEdge(elem2edge(bdFlag(:,1) == 1,[4,5,6])) = true;
        isFixedEdge(elem2edge(bdFlag(:,2) == 1,[2,3,6])) = true;
        isFixedEdge(elem2edge(bdFlag(:,3) == 1,[1,3,5])) = true;
        isFixedEdge(elem2edge(bdFlag(:,4) == 1,[1,2,4])) = true;
        bdEdge = edge(isFixedEdge,:);
        isFixedNode = false(N,1);
        isFixedNode(bdEdge(:)) = true;
    end
    fixedNode = find(isFixedNode);
    fixedEdge = find(isFixedEdge);
    fixedDof = [fixedNode; N + fixedEdge];
    freeNode = find(~isFixedNode);
    freeEdge = find(~isFixedEdge);
    freeDof = [freeNode; N + freeEdge];
    isPureNeumann = false;
    if isempty(fixedNode) % pure Neumann boundary condition
        isPureNeumann = true;
        fixedDof = Ndof;
        freeDof = (1:Ndof-1)';    % eliminate the kernel
    end

    % nonzero Dirichlet boundary condition
    if ~isfield(pde,'gu') || (isnumeric(pde.gu) && all(pde.gu == 0))   % zero gu
        pde.gu = [];
    end
    if ~isfield(pde,'gsigma') || (isnumeric(pde.gsigma) && all(pde.gsigma == 0))   % zero gsigma
        pde.gsigma = [];
    end
    if ~isPureNeumann && ~isempty(fixedNode) && (~isempty(pde.gu) || ~isempty(pde.gsigma))
        if (isnumeric(pde.gu) && length(pde.gu) == NE)
            u(isFixedEdge) = pde.gu(isFixedEdge);
        else
            u(isFixedEdge) = edgeinterpolate(pde.gu,node,bdEdge);
        end
        sigma(fixedNode) = pde.gsigma(node(fixedNode,:));
        f = f - A*[sigma; u];
        f(fixedDof) = [sigma(fixedNode); u(fixedEdge)];    
    end
    % modify the matrix to include the Dirichlet boundary condition
    if ~isempty(fixedDof)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedDof) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end 
    
    % Pure Neumann boundary condition
    if isPureNeumann
        f(N+1:Ndof) = f(N+1:Ndof) - mean(f(N+1:Ndof));   % compatilbe condition: sum(b) = 0
    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions bdfaceintegral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function bdb = bdfaceintegral(isBdElem,face,face2locdof)
    %% Compute boundary surface integral of lowest order edge element.
    %  bdb(k) = \int_{face} (n×gcurlu, phi_k) dS

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
        gNp = pde.gcurlu(pxyz,normal);    
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
            idx(:,s) = elem2edge(isBdElem,kk);
        end
    end
    %% Distribute to DOF
    bdb = accumarray(idx(:),bt(:),[Ndof 1]);        
    end

end
