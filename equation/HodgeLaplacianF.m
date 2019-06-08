function [sigma,u,AD] = HodgeLaplacianF(node,elem,bdFlag,pde,option)
%% HODEGELAPLACIANF Hodge Laplacian of edge element


if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

%% Data structure
% elemold = elem;
[elem,bdFlag] = sortelem(elem,bdFlag);  % ascend ordering
[elem2edge,edge] = dofedge(elem);
[Dlambda,area,elemSign] = gradbasis(node,elem);
locEdge = [2 3; 1 3; 1 2];
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
Nsigma = N; Nu = NE; Ndof = Nsigma + Nu;

%% Assemble matrix 
% Mass matrices
Mv = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
Mv = spdiags(Mv,0,N,N);
Me = getmassmatvec(elem2edge,area,Dlambda,'RT0');
Mt = spdiags(1./area,0,NT,NT);
% C. curl operator: P1 -> (RT0)'
C = Me*icdmat(double(edge),[-1 1]);  % gradient matrix
% R. rotation operator
B = icdmat(double(elem2edge),[1 -1 1]);
D = B'*Mt*B;
A = [-Mv C'; C D];

%% Assemble right hand side
f = zeros(Ndof,1);
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
            % phi_k = rotation(lambda_iDlambda_j - lambda_jDlambda_i);
            phi_k(:,1) = -(lambda(p,i)*Dlambda(:,2,j)-lambda(p,j)*Dlambda(:,2,i));
            phi_k(:,2) = lambda(p,i)*Dlambda(:,1,j)-lambda(p,j)*Dlambda(:,1,i);
            rhs = dot(phi_k,fp,2);
            bt(:,k) = bt(:,k) + w(p)*rhs;
        end
    end
    bt = bt.*repmat(area,1,3);
    f = accumarray(N+elem2edge(:),bt(:),[Ndof 1]);
end
clear pxy fp bt rhs phi_k Dlambda

%% Boundary Conditions
[AD,f,u,sigma,freeDof,isPureNeumann] = getbdHodgeLapE(f);

%% Solve the linear system
temp = zeros(Ndof,1);
% temp(freeDof) = A(freeDof,freeDof)\f(freeDof);
AD = A(freeDof,freeDof);
sigma(freeNode) = temp(freeNode);
u(freeEdge) = temp(freeEdge+N);

if isPureNeumann
   u = u - mean(u); % normalize u 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions getbdHodgeLapE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,f,u,sigma,freeDof,isPureNeumann] = getbdHodgeLapE(f)

    u = zeros(NE,1); sigma = zeros(N,1);    
    %% Neumann boundary condition
    % Find boundary edges: Neumann
    isNeumann(elem2edge(bdFlag(:)==2)) = true;
    Neumann = edge(isNeumann,:); 
    % Part sigma: int_e (gut, tau)
    if isnumeric(pde.gut) && all(pde.gut == 0)
        pde.gut = [];
    end
    if ~isempty(Neumann) && ~isempty(pde.gut)
        el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [lambdagN,weightgN] = quadpts1(option.gNquadorder);
        phigN = lambdagN;                 % linear bases
        nQuadgN = size(lambdagN,1);
        ge = zeros(size(Neumann,1),2);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(pp,2)*node(Neumann(:,2),:);
            gNp = pde.gut(ppxy);
            for igN = 1:2
                ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
            end
        end
        ge = ge.*repmat(el,1,2);
        f = f + accumarray(Neumann(:), ge(:),[Ndof,1]); 
    end
    % Part u: int_F (div u, v dot n)
    edgeSign = ones(NE,1);
    idx = (bdFlag(:,1) ~= 0) & (elemSign == -1);% first edge is on boundary
    edgeSign(elem2edge(idx,1)) = -1;
    idx = (bdFlag(:,2) ~= 0) & (elemSign == 1); % second edge is on boundary
    edgeSign(elem2edge(idx,2)) = -1;
    idx = (bdFlag(:,3) ~= 0) & (elemSign == -1);% first edge is on boundary
    edgeSign(elem2edge(idx,3)) = -1;
    if isnumeric(pde.gdivu) && all(pde.gdivu==0)
        pde.gdivu = [];
    end
    if ~isempty(pde.gdivu) && ~isempty(Neumann)
        if ~isfield(option,'gNquadorder')
            option.gNquadorder = 2;   % default order exact for linear gN
        end
        [lambda,weight] = quadpts1(option.gNquadorder);
        nQuad = size(lambda,1);
        Neumannidx = find(isNeumann) + N;
        for ip = 1:nQuad
        	pxy = lambda(ip,1)*node(Neumann(:,1),:)+...
                  lambda(ip,2)*node(Neumann(:,2),:);               
            f(Neumannidx) = f(Neumannidx) + weight(ip)*pde.gdivu(pxy);
        end
        f(Neumannidx) = f(Neumannidx).*edgeSign(isNeumann);
        % no edge length since the basis of edge element cancle it.
    end
    
    % Find Dirichlet boundary nodes and edges
    isDirichlet = false(NE,1);
    isDirichlet(elem2edge(bdFlag(:)==1)) = true;
    Dirichlet = edge(isDirichlet,:);
    isfixedNode = false(N,1); 
    isfixedNode(Dirichlet(:)) = true;
    fixedNode = find(isfixedNode);
    fixedDof = [fixedNode; N + find(isDirichlet)];
    freeNode = find(~isfixedNode);
    freeEdge = find(~isDirichlet);
    freeDof = [freeNode; N + freeEdge];
    isPureNeumann = false;
    if isempty(fixedNode) % pure Neumann boundary condition
        isPureNeumann = true;
        fixedDof = Ndof;
        freeDof = (1:Ndof-1)';    % eliminate the kernel
    end
    
    % Modify right hand side to include Dirichlet boundary condition
    if isnumeric(pde.gu) && all(pde.gu == 0)   % zero gu
        pde.gu = [];
    end
    if isnumeric(pde.gsigma) && all(pde.gsigma == 0)   % zero gsigma
        pde.gsigma = [];
    end
    if ~isPureNeumann && ~isempty(fixedNode) && (~isempty(pde.gu) || ~isempty(pde.gsigma))
        sigma = zeros(N,1);
        sigma(fixedNode) = pde.gsigma(node(fixedNode,:));
        u = zeros(NE,1);
        u(isDirichlet) = edgeinterpolate(pde.gu,node,Dirichlet);
        f = f - A*[sigma; u];
        f(fixedDof) = [sigma(fixedNode); u(isDirichlet)];    
    end
    
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedNode,fixedNode)=I, AD(fixedNode,freeNode)=0, AD(freeNode,fixedNode)=0.
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
end
