function [u,w,AE,AI,eqn,info] = interfacePoisson(node,elem,pde,interfaceEdge,option)
%% INTERFACEPOISSON Poisson equation: P1 linear element.
%
%   u = INTERFACEPOISSON(node,elem,pde,E) produces the linear finite element
%   approximation of the interface Poisson equation
% 
%       -div(d*grad(u))=f  in \Omega, with 
%       Dirichlet boundary condition u=g_D on \Gamma_D, 
%       Neumann boundary condition   d*grad(u)*n=g_N on \Gamma_N,
%       Robin boundary condition     g_R*u + d*grad(u)*n=g_N on \Gamma _R
%
% Huayi: type a brief describtion of the problem.

N = size(node,1);	NT = size(elem,1);    Ndof = N;

%% Geometry structures of interface meshes

isExteriorElem = false(NT,1);
center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
isExteriorElem(pde.phi(center)>0) = true;

%% Diffusion coefficient
if isfield(pde,'d') && ~isempty(pde.d)
   if isnumeric(pde.d)
      K = pde.d;                   % d is an array
   else                            % d is a function
      K = pde.d(center);  
   end
else
    K = [];
end

%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
A = sparse(N,N);
AE = sparse(N,N);
AI = sparse(N,N);
for i = 1:3
    for j = i:3
        Aij = (Dphi(:,1,i).*Dphi(:,1,j)+Dphi(:,2,i).*Dphi(:,2,j)).*area;
        if ~isempty(K)
            Aij = K.*Aij;
        end
        if (j==i)
            A = A + sparse(elem(:,i),elem(:,j),Aij,N,N);
            AE = AE + sparse(elem(isExteriorElem,i),elem(isExteriorElem,j),Aij(isExteriorElem),N,N);
            AI = AI + sparse(elem(~isExteriorElem,i),elem(~isExteriorElem,j),Aij(~isExteriorElem),N,N);
        else
            A = A + sparse([elem(:,i);elem(:,j)],[elem(:,j);elem(:,i)],[Aij; Aij],N,N);  
            AE = AE + sparse([elem(isExteriorElem,i);elem(isExteriorElem,j)],...
                [elem(isExteriorElem,j);elem(isExteriorElem,i)],[Aij(isExteriorElem);Aij(isExteriorElem)], N,N);
            AI = AI + sparse([elem(~isExteriorElem,i);elem(~isExteriorElem,j)],...
                [elem(~isExteriorElem,j);elem(~isExteriorElem,i)],[Aij(~isExteriorElem);Aij(~isExteriorElem)], N,N);
        end        
    end
end
clear  Aij

%% Assemble the right hand side
b = zeros(Ndof,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    [lambda,weight] = quadpts(option.fquadorder);
    nQuad = size(lambda,1);
    ft = zeros(NT,3);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        % function values at quadrature points
        fp = pde.f(pxy);
        % evaluate fp outside.
        for j = 1:3
            ft(:,j) = ft(:,j) + lambda(p,j)*weight(p)*fp;
        end
    end
    ft = ft.*[area,area,area];
    b = accumarray(elem(:),ft(:),[Ndof 1]);
end

%% Extend w to the whole domain
isInterfaceNode = false(N,1);
isInterfaceNode(interfaceEdge(:)) = true;
interfaceNode = find(isInterfaceNode);

isInNode = false(N,1);
interiorElem = elem(~isExteriorElem,:);
isInNode(interiorElem(:)) = true;
inNode = find(isInNode & ~isInterfaceNode);

w = zeros(N,1);
FI= zeros(N,1);
w(interfaceNode) = pde.exactw(node(interfaceNode,:));
FI = FI - AI*w;
extensionoption.printlevel = 0;
w(inNode) = amg(AI(inNode, inNode),FI(inNode),extensionoption);

%% Neumann boundary condition on interface edges
ve = node(interfaceEdge(:,1),:) - node(interfaceEdge(:,2),:);
ve = [-ve(:,2), ve(:,1)];
edgeLen = sqrt(sum(ve.^2,2));
if ~isfield(option,'gNquadorder')
    option.gNquadorder = 5;   % default order
end
[lambda,weight] = quadpts1(option.gNquadorder);
nQuad = length(weight);
ge = zeros(size(interfaceEdge,1),2);
for i = 1:nQuad
    pxy=lambda(i,1)*node(interfaceEdge(:,1),:)+lambda(i,2)*node(interfaceEdge(:,2),:);
    leftpt = pxy - ve/2;
    rightpt = pxy + ve/2;
    pxy = findintersectbisect(pde.phi,leftpt,rightpt);
    q = pde.exactq(pxy);
    ge(:,1)=ge(:,1)+weight(i)*lambda(i,1)*q;
    ge(:,2)=ge(:,2)+weight(i)*lambda(i,2)*q;
end
ge = ge.*[edgeLen,edgeLen];
b = b - accumarray(interfaceEdge(:), ge(:),[Ndof,1]);

%% Dirichlet boundary condition
[~,bdEdge] = findboundary(elem);
isBdNode = false(N,1); 
isBdNode(bdEdge(:)) = true;
bdNode = find(isBdNode);

u = zeros(N,1); 
u(bdNode) = pde.g_D(node(bdNode,:));
b = b - A*u+AI*w;
b(bdNode) = u(bdNode);

freeNode = find(~isBdNode);
bdidx = zeros(N,1); 
bdidx(bdNode) = 1;
Tbd = spdiags(bdidx,0,N,N);
T = spdiags(1-bdidx,0,N,N);
AD = T*A*T + Tbd;

%% Solve the system of linear equations
if isempty(freeNode), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 1e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option);                 
end

%%
eqn = struct('A',AD,'b',b,'freeNode',freeNode);