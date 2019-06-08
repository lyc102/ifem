function [sigma,u] = HodgeLaplacian3F(node,elem,bdFlag,pde,option)
%% HODEGELAPLACIAN3E Hodge Laplacian of face element


if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

%% Data structure
% elemold = elem;
N = size(node,1); NT = size(elem,1); 
[elem,bdFlag] = sortelem3(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dof3edge(elem);
[elem2face,face] = dof3face(elem);
NF = size(face,1); NE = size(edge,1);
Nsigma = NE; Nu = NF; Ndof = Nsigma + Nu;
face2edge = zeros(NF,3,'int32');
face2edge(elem2face(:,1),:) = elem2edge(:,[4 5 6]);
face2edge(elem2face(:,2),:) = elem2edge(:,[2 3 6]);
face2edge(elem2face(:,3),:) = elem2edge(:,[1 3 5]);
face2edge(elem2face(:,4),:) = elem2edge(:,[1 2 4]);

[Dlambda,volume,elemSign] = gradbasis3(node,elem);
localFace = [2 3 4; 1 3 4; 1 2 4; 1 2 3]; 

%% Assemble matrix 
% Mass matrices
Me = getmassmatvec3(elem2edge,volume,Dlambda,'ND0');
Mf = getmassmatvec3(elem2face,volume,Dlambda,'RT0');
Mt = spdiags(1./volume,0,NT,NT);
% C. curl operator: ND0 -> (RT0)'
C = Mf*icdmat(double(face2edge),[1 -1 1]);  % curl matrix
% B'MB: -grad div operator
B = icdmat(double(elem2face),[1 -1 1 -1]);  % divergence matrix
D = B'*Mt*B;
% big matrix for the mixed system
A = [-Me C'; C D];

%% Assemble right hand side
f = zeros(Ndof,1);
if ~isfield(pde,'f') || (isfield(pde,'f') && isreal(pde.f) && all(pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f)
    fI = faceinterpolate3(pde.f,node,face,2);
    f = [zeros(NE,1); Mf*fI];
end

%% Boundary Conditions
[AD,f,u,sigma,freeDof,freeEdge,freeFace,isPureNeumann] = getbdHodgeLap3F(f);

%% Solve the linear system
temp = zeros(Ndof,1);
temp(freeDof) = A(freeDof,freeDof)\f(freeDof);
sigma(freeEdge) = temp(freeEdge);
u(freeFace) = temp(freeFace+NE);

if isPureNeumann
   u = u - mean(u); % normalize u 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdHodgeLapE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,f,u,sigma,freeDof,freeEdge,freeFace,isPureNeumann] = getbdHodgeLap3F(f)

    %% Neumann boundary condition
    isNeumann(elem2face(bdFlag(:)==2)) = true;
    Neumann = face(isNeumann,:); 
    % Part sigma: int_F (u cross n, tau)
    g = zeros(Ndof,1);
    if isnumeric(pde.gu) && all(pde.gu == 0)   % zero gu
        pde.gu = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.gu)
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
    % Part u: int_F (div u, v dot n)
    faceSign = ones(NF,1);
    % elementwise consistence sign:    [1 -1 1 -1]
    idx = (bdFlag(:,1) ~= 0) & (elemSign == -1);% first face is on boundary
    faceSign(elem2face(idx,1)) = -1;
    idx = (bdFlag(:,2) ~= 0) & (elemSign == 1); % second face is on boundary
    faceSign(elem2face(idx,2)) = -1;
    idx = (bdFlag(:,3) ~= 0) & (elemSign == -1);% third face is on boundary
    faceSign(elem2face(idx,3)) = -1;    
    idx = (bdFlag(:,4) ~= 0) & (elemSign == 1);% third face is on boundary
    faceSign(elem2face(idx,4)) = -1;    
    if isnumeric(pde.gdivu) && all(pde.gdivu == 0)
        pde.gdivu = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.gdivu)
        v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
        v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        % three middle points rule
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
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
            f(Neumannidx) = f(Neumannidx) + weight(ip)*pde.gdivu(ppxyz);
        end
        f(Neumannidx) = f(Neumannidx).*faceSign(isNeumann);
        % no area length since the basis of face element cancle it.        
    end

    %% Dirichlet boundary condition
    % Find Dirichlet boundary dof: fixedDof
    isDirichlet = false(NF,1);
    isDirichlet(elem2face(bdFlag(:)==1)) = true;
    Dirichlet = face(isDirichlet,:);
    fixedFace = find(isDirichlet);
    isFixedEdge = false(NE,1);
    isFixedEdge(face2edge(isDirichlet,:)) = true;
    fixedEdge = find(isFixedEdge);
    fixedDof = [fixedEdge; NE + fixedFace];
    freeEdge = find(~isFixedEdge);
    freeFace = find(~isDirichlet);
    freeDof = [freeEdge; NE + freeFace];
    isPureNeumann = false;
    if isempty(fixedFace) % pure Neumann boundary condition
        isPureNeumann = true;
        fixedDof = Ndof;
        freeDof = (1:Ndof-1)';    % eliminate the kernel
    end

    % nonzero Dirichlet boundary condition
    if isnumeric(pde.gu) && all(pde.gu == 0)   % zero gu
        pde.gu = [];
    end
    u = zeros(NF,1);
    if ~isPureNeumann && ~isempty(fixedFace) && ~isempty(pde.gu)
        if (isnumeric(pde.gu) && length(pde.gu) == NF)
            u(isDirichlet) = pde.gu(isDirichlet);
        else
            u(isDirichlet) = faceinterpolate3(pde.gu,node,Dirichlet);
        end
    end
    if isnumeric(pde.gsigma) && all(pde.gsigma == 0)   % zero gsigma
        pde.gsigma = [];
    end
    sigma = zeros(NE,1);
    if ~isPureNeumann && ~isempty(fixedFace) && ~isempty(pde.gsigma)
        bdEdge = edge(isFixedEdge,:);
        sigma(isFixedEdge) = edgeinterpolate(pde.gsigma,node,bdEdge);
    end
    f = f - A*[sigma; u];
    f(fixedDof) = [sigma(fixedEdge); u(fixedFace)];    

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
    %% Compute boundary surface integral of lowest order face element.
    %  bdb(k) = \int_{face} (n×gcurlu, phi_k) dS

    %% Compute scaled normal
    faceIdx = true(4,1);
    faceIdx(face) = false;
    normal = -3*repmat(volume(isBdElem),1,3).*Dlambda(isBdElem,:,faceIdx);

    %% Data structure
    tetLocEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % face of a tetrahedral [1 2 3 4]
    face2locEdge = [2 3; 3 1; 1 2]; % face of the face [1 2 3]

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
        gNp = mycross(pde.gu(pxyz),normal,2);    
        for s = 1:3
            kk = face2locdof(s);
            pidx = face(face2locEdge(s,1))< face(face2locEdge(s,2));
            % phi_k = lambda_iDlambda_j - lambda_jDlambda_i;
            % lambda_i is associated to the local index of the face [1 2 3]
            % Dlambda_j is associtated to the index of tetrahedron
            % - when the direction of the local face s is consistent with the
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
            bt(:,s) = bt(:,s) + w(pp)*rhs; % area is included in normal
            idx(:,s) = elem2face(isBdElem,kk);
        end
    end
    %% Distribute to DOF
    bdb = accumarray(idx(:),bt(:),[Ndof 1]);        
    end

end
