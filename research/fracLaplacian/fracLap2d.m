function [u,eqn,info,node,elem,Neumann] = fracLap2d(cube,h,pde,option)

if ~exist('option','var'), option = []; end
%% Parameters
global s
alpha = 1-2*s;
if s == 0.5
    gamma = 1;
else
    gamma = 3/(2*s)+0.1;
end

if ~exist('option','var'), option = []; end
if ~isfield(option,'plotflag'), option.plotflag = 0; end

%% Mesh
% we still use y to denote the one additional dimension
if option.gradmesh
   [node,elem,y] = cubehexgradmesh(cube,h,gamma,option.plotflag);
else
    [node,elem] = cubehexmesh(cube,h);
    y0 = cube(5); y1 = cube(6);
    y = (y0:h:y1)';
end
N = size(node,1); NT = size(elem,1); Ndof = N;

%% Stiffness matrix and mass matrix in each direction
tic;
% y direction
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
% x-direction
hx1 = h; % node(elem(:,2),1) - node(elem(:,1),1);
hx2 = h; % node(elem(:,4),2) - node(elem(:,1),2);
Ax1 = 1./hx1*[1 -1];
Mx1 = hx1*[1/3 1/6];
Ax2 = 1./hx2*[1 -1];
Mx2 = hx2*[1/3 1/6];

%% Assemble stiffness matrix
% generate sparse pattern
ii = zeros(36*NT,1); jj = zeros(36*NT,1); 
index = 0;
for i = 1:8
    for j = i:8
        ii(index+1:index+NT) = double(elem(:,i)); 
        jj(index+1:index+NT) = double(elem(:,j));  
        index = index + NT;
    end
end
% index map
itoixiyiz = [[1 2 2 1 1 2 2 1]' [1 1 2 2 1 1 2 2]' [1 1 1 1 2 2 2 2]'];
ij2k(1,1) = 1; ij2k(1,2) = 2; ij2k(2,1) = 2; ij2k(2,2) = 3; % index map for Ay 

if ~isfield(pde,'d'), pde.d = []; end
sA = zeros(36*NT,1);
index = 0;
for i = 1:8
    for j = i:8
        ix1 = itoixiyiz(i,1); jx1 = itoixiyiz(j,1);
        ix2 = itoixiyiz(i,2); jx2 = itoixiyiz(j,2);
        iy = itoixiyiz(i,3); jy = itoixiyiz(j,3);
        Aij = Ax1(2-(ix1==jx1))*Mx2(2-(ix2==jx2)).*My(:,ij2k(iy,jy)) ...
            + Mx1(2-(ix1==jx1))*Ax2(2-(ix2==jx2)).*My(:,ij2k(iy,jy)) ...
            + Mx1(2-(ix1==jx1))*Mx2(2-(ix2==jx2)).*Ay(:,ij2k(iy,jy));
        sA(index+1:index+NT) = repmat(Aij,NTx,1);
        index = index + NT;
    end
end
% assemble the matrix
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj sA

%% Assemble the right hand side
% harmonic extension
b = zeros(Ndof,1);
u = zeros(Ndof,1);

%% Set up boundary conditions
nx1 = (cube(2)-cube(1))/h+1;  % number of grid points in x-direction
nx2 = (cube(4)-cube(3))/h+1;  % number of grid points in x-direction
ny = length(y);  % number of grid points in y-direction
nodeidx = reshape(1:N,ny,nx1,nx2);
isBdNode = true(N,1);
freeNode = nodeidx(1:ny-1,2:nx1-1,2:nx2-1);
freeNode = freeNode(:);
isBdNode(freeNode(:)) = false;

% Modify the matrix to include the Dirichlet boundary condition
bdidx = zeros(Ndof,1); 
bdidx(isBdNode) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
A = T*A*T + Tbd;

% Neumann faces
bottom = squeeze(nodeidx(1,:,:));
nodeidx2d = reshape(1:length(bottom(:)),nx1,nx2);
k = nodeidx2d(1:nx1-1,1:nx2-1);
k = k(:);
% couple with the index map. bottom is the orignal index
Neumann = [bottom(k) bottom(k+1) bottom(k+nx1+1) bottom(k+nx1)];
area = hx1.*hx2;

% compute integral
if ~isfield(option,'gNquadorder')
    option.gNquadorder = 3;   % default order exact for linear gN
end
[GaussPt,weightgN] = quadptsquad(option.gNquadorder);
lambdax1 = GaussPt(:,1);
lambdax2 = GaussPt(:,2);
nQuadgN = size(GaussPt,1);
% barycentrical coordinate for four basis
ge = zeros(size(Neumann,1),4);
for pp = 1:nQuadgN
    phi = kron([lambdax2(pp) 1-lambdax2(pp)],[lambdax1(pp); 1-lambdax1(pp)]);
    phi = phi([1 2 4 3]); % switch the index of 3 and 4
    % quadrature points in the x-y coordinate
    ppxyz = node(Neumann(:,1),1:2)*phi(1) + node(Neumann(:,2),1:2)*phi(2)...
          + node(Neumann(:,3),1:2)*phi(3) + node(Neumann(:,4),1:2)*phi(4);
    gNp = pde.g_N(ppxyz);
    for iN = 1:4
        ge(:,iN) = ge(:,iN) + weightgN(pp)*phi(iN)*gNp;
    end
end
ge = ge*area;
b = b + accumarray(Neumann(:),ge(:),[Ndof,1]); 
% Neumann edges are considered as open set. So the corner points should be
% set as Dirichlet boundary condition!
b(isBdNode) = 0;

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeNode), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
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
        u(freeNode) = A(freeNode,freeNode)\b(freeNode);
        residual = norm(b - A*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
%         option.x0 = u;
        option.solver = 'CG';
%         [u,info] = mg(A,b,elem,option);
        [u(freeNode),info] = amg(A(freeNode,freeNode),b(freeNode),option);                 
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(A(freeNode,freeNode),b(freeNode),option);                 
end

%% Compute error using boundary integral
errH1 = getH1error3bd(node,Neumann,pde,u,A,5);

%% Output information
eqn = struct('A',A,'b',b,'freeNode',freeNode);
info.errH1 = errH1;
info.assembleTime = assembleTime;