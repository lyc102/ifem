function [u,eqn,info,node,elem,Neumann] = fracLap1d(square,h,pde,option)

global s
alpha = 1-2*s;
if s == 0.5
    gamma = 1;
else
    gamma = 3/(2*s)+0.1;
end

tic;
%% Mesh
if option.gradmesh
   [node,elem,y] = squaregradmeshquad(square,h,gamma);   
else
   [node,elem] = squarequadmesh(square,h);
    y0 = square(3); y1 = square(4);
    y = (y0:h:y1)';
end
N = size(node,1); NT = size(elem,1); Ndof = N;
% x0 = square(1); x1 = square(2); 

%% Compute Local matrix
hy = diff(y);
NTy = length(hy);
NTx = NT/NTy;
% stiffness matrix in y direction
Ay = zeros(NTy,3);
intyalpha = diff(y.^(alpha+1)/(alpha+1));
Ay(:,1) = intyalpha./hy.^2;  % Ay(1,1)
Ay(:,2) = -Ay(:,1);          % Ay(1,2) and Ay(2,1)
Ay(:,3) = Ay(:,1);           % Ay(2,2)
% mass matrix in y direction
My = zeros(NTy,3);
a3 = diff(y.^(alpha+3)/(alpha+3));
a2 = diff(y.^(alpha+2)/(alpha+2));
a1 = diff(y.^(alpha+1)/(alpha+1));
yleft = y(1:end-1);
yright = y(2:end);
My(:,1) = (a3 - 2*yleft.*a2 + yleft.^2.*a1)./hy.^2;  % My(1,1)
My(:,2) = (-a3 + (yleft+yright).*a2 - yleft.*yright.*a1)./hy.^2;    % My(1,2) and My(2,1)
My(:,3) = (a3 - 2*yright.*a2 + yright.^2.*a1)./hy.^2;  % My(2,2)
% stiffness matrix in x direction
hx = h;
Ax = 1/hx*[1 -1; -1 1];
Mx = hx/6*[2 1; 1 2];

%% Assemble stiffness matrix
% index map 
ixiy = [1 1; 2 1; 2 2; 1 2];
ij2k(1,1) = 1; ij2k(1,2) = 2; ij2k(2,1) = 2; ij2k(2,2) = 3;
% generate sparse pattern
ii = zeros(10*NT,1); jj = zeros(10*NT,1); 
index = 0;
for i = 1:4
    for j = i:4
        ii(index+1:index+NT) = double(elem(:,i)); 
        jj(index+1:index+NT) = double(elem(:,j));  
        index = index + NT;
    end
end
% compute non-zeros
sA = zeros(10*NT,1);
index = 0;
for i = 1:4
    for j = i:4
        ix = ixiy(i,1); iy = ixiy(i,2);
        jx = ixiy(j,1); jy = ixiy(j,2);        
        Aij = Ax(ix,jx)*My(:,ij2k(iy,jy)) + Mx(ix,jx)*Ay(:,ij2k(iy,jy));
        sA(index+1:index+NT,1) = repmat(Aij,NTx,1);
        index = index + NT;
    end
end    

% assemble the matrix
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj AU

%% Right hand side 
% f = 0
b = zeros(Ndof,1);

%% Set up boundary conditions
% find boundary nodes and Neumann edges
nx = NTx + 1;
ny = NTy + 1;
fixedNode = [1:ny (2:nx-1)*ny (nx-1)*ny+(1:ny)];
bottom = transpose(1:ny:(nx-1)*ny+1);
Neumann = [bottom(1:end-1) bottom(2:end)];
isBdNode = false(N,1);
isBdNode(fixedNode) = true;
freeNode = find(~isBdNode);

% 
% bdFlag = setboundary(node,elem,'Dirichlet','abs(y)>eps','Neumann','abs(y)<=eps');
% [fixedNode,Neumann,isBdNode] = findboundary(elem,bdFlag);

% Modify the matrix to include the Dirichlet boundary condition
bdidx = zeros(Ndof,1); 
bdidx(fixedNode) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
A = T*A*T + Tbd;

% Neumann boundary condition
el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
if ~isfield(option,'gNquadorder')
    option.gNquadorder = 4;   % default order exact for linear gN
end
[lambdagN,weightgN] = quadpts1(option.gNquadorder);
phigN = lambdagN;                 % linear bases
nQuadgN = size(lambdagN,1);
ge = zeros(size(Neumann,1),2);
for pp = 1:nQuadgN
    % quadrature points in the x-y coordinate
    ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
         + lambdagN(pp,2)*node(Neumann(:,2),:);
    gNp = pde.g_N(ppxy);
    for igN = 1:2
        ge(:,igN) = ge(:,igN) + weightgN(pp)*phigN(pp,igN)*gNp;
    end
end
ge = ge.*repmat(el,1,2);
b = b + accumarray(Neumann(:), ge(:),[Ndof,1]); 

% Modify right handside for Dirichlet boundary condition
% Neumann edges are considered as open set. So the corner points should be
% set as Dirichlet boundary condition!
b(fixedNode) = 0;

eqn = struct('A',A,'b',b,'freeNode',freeNode);

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve
u = zeros(Ndof,1);
switch option.solver
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
        return
    case 'direct'
        u(freeNode) = A(freeNode,freeNode)\b(freeNode);
        residual = norm(b - A*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);        
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(A(freeNode,freeNode),b(freeNode),option);                 
end
info.assembleTime = assembleTime;

%% Compute error using boundary integral
[lambdagN,weightgN] = quadpts1(option.gNquadorder);
phigN = lambdagN;                 % linear bases
nQuadgN = size(lambdagN,1);
err = zeros(size(Neumann,1),1);
for pp = 1:nQuadgN
    % quadrature points in the x-y coordinate
    ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
         + lambdagN(pp,2)*node(Neumann(:,2),:);
    uhp = u(Neumann(:,1))*phigN(pp,1) + u(Neumann(:,2))*phigN(pp,2);
    up = pde.exactu(ppxy);
    gNp = pde.g_N(ppxy);
    err = err + weightgN(pp)*gNp.*(up - 2*uhp);
end
err = sum(err.*el) + u'*A*u;
errH1 = sqrt(abs(err));
info.errH1 = errH1;