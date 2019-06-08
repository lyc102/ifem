function [u,eqn,info] = fracLapP1P1(node,elem,pde,option)
%% FRACLAPP1P1 solves fractional Laplacian equation using P1-P1 element
%
% [u,eqn,info] = fracLapP1P1(node,elem,pde,option) solves the fractional
% Laplacian equation
% 
%  (-\Delta^^s u = f in \Omega with u = 0 on \partial \Omega.
%
% 
%
%% Options
if ~exist('option','var'), option = []; end
if ~isfield(option,'plotflag'), option.plotflag = 0; end
if ~isfield(option,'yrefinement'), option.yrefinement = 0; end

%% Parameters
global s
alpha = 1-2*s;
if s == 0.5
    Gamma = 1;
else
    Gamma = 3/(2*s)+0.01;
%     Gamma = 7/(2*s)+0.01;
end

%% Mesh
N = size(node,1); 
NT = size(elem,1);
L = 1+log10(NT)/3; % to control the truncated error, L should increase as logN
% L = pde.L;
My = round(L*sqrt(NT)); % number of elements in y direction
y = gradmap(0,L,Gamma,My);
if option.yrefinement % further refinement to enforce mesh condition (5.1)
    area = simplexvolume(node,elem);
    hxmin = min(sqrt(area));
    hy = diff(y);
    idx = find(hy > option.yrefinement*hxmin);
    while any(idx)
        ynew = (y(idx) + y(idx+1))/2;
        y = [y; ynew]; %#ok<AGROW>
        y = sort(y);
        hy = diff(y);
        idx = find(hy > option.yrefinement*hxmin);
    end
end
Ny = length(y); % number of vertices in y-direction  
NTy = Ny - 1;   % number of elements in y-direction
NTtotal = NT*NTy;% total number of elements of the tensor product mesh
Ntotal = N*Ny;   % total number of vertices of the tensor product mesh
Ndof = Ntotal;   % total degree of freedoms

tic; % start the timing for the assembling
%% Stiffness matrix and mass matrix in the extended direction
% quantities in y direction
hy = diff(y);
a = zeros(NTy,3);
for i = 1:3
    a(:,i) = diff(y.^(alpha+i)/(alpha+i));
end
% stiffness matrix in y direction
Ay = zeros(NTy,2,2);
Ay(:,1,1) = a(:,1)./hy.^2;
Ay(:,1,2) = -Ay(:,1,1);           
Ay(:,2,1) = -Ay(:,1,1);           
Ay(:,2,2) = Ay(:,1,1);            
% mass matrix in y direction
My = zeros(NTy,2,2);
y1 = y(1:end-1);
y2 = y(2:end);
My(:,1,1) = (a(:,3) - 2*y1.*a(:,2) + y1.^2.*a(:,1))./hy.^2;
My(:,1,2) = (-a(:,3) + (y1+y2).*a(:,2) - y1.*y2.*a(:,1))./hy.^2;
My(:,2,1) = My(:,1,2);
My(:,2,2) = (a(:,3) - 2*y2.*a(:,2) + y2.^2.*a(:,1))./hy.^2;

%% Stiffness matrix and mass matrix in the original direction
[Dphi,area] = gradbasis(node,elem); 
% Compute a piecewise constant diffusion coefficient
if ~isfield(pde,'d'), pde.d = []; end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end
if ~isempty(pde.d) && isnumeric(pde.d)
   K = pde.d;                                 % d is an array
end
if ~isempty(pde.d) && ~isnumeric(pde.d)       % d is a function   
    [lambda,weight] = quadpts(option.dquadorder);
    nQuad = size(lambda,1);
    K = zeros(NT,1);
    for p = 1:nQuad
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
        K = K + weight(p)*pde.d(pxy);      
   end
end
if ~isempty(pde.d) % build the coefficients into the scaled area
    areaK = K.*area;
else
    areaK = area;
end
At = zeros(NT,3,3);
Mt = zeros(NT,3,3);
for i = 1:3
    for j = 1:3
        At(:,i,j) = (Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j)).*areaK;
        Mt(:,i,j) = areaK*((i==j)+1)/12;
    end
end

%% Assemble the stiffness matrix
dofMap = repmat(1:N,Ny,1) + repmat((0:NTy)'*N,1,N);
ii = zeros(36*NTtotal,1); 
jj = zeros(36*NTtotal,1); 
sA = zeros(36*NTtotal,1);
index = 0;
for m = 1:2
    for n = 1:2
        for i = 1:3
            for j = 1:3
                ii(index+1:index+NTtotal) = dofMap(m:Ny+m-2,elem(:,i));
                jj(index+1:index+NTtotal) = dofMap(n:Ny+n-2,elem(:,j));
                sA(index+1:index+NTtotal) = kron(At(:,i,j),My(:,m,n)) + ...
                                            kron(Mt(:,i,j),Ay(:,m,n));
                index = index + NTtotal;
            end
        end
    end
end
A = sparse(ii,jj,sA,Ndof,Ndof);
clear ii jj sA At Mt Ay My

%% Compute average aspect ratio
hx = mean(sqrt(area));
% aspectRatio = hx./hy(end:-1:1);
aspectRatio = hx/hy(1);
% hxinv = mean(1./sqrt(area));
% aspectRatio = hy(end:-1:1)*hxinv;

%% Assemble the right hand side
b = zeros(Ndof,1);
u = zeros(Ndof,1);

%% Set up boundary conditions
% find Dirichlet boundary nodes
bdNode = findboundary(elem);  % boundary vertices of mesh in x-direction
lateralbd = dofMap(1:Ny-1,bdNode);
topbd = dofMap(Ny,:);

% Modify the matrix to include the Dirichlet boundary condition
bdidx = zeros(Ndof,1); 
bdidx(lateralbd) = 1;
bdidx(topbd) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
A = T*A*T + Tbd;

% Compute boundary integral over the original domain
if ~isfield(option,'fquadorder')
    option.fquadorder = 4;   
end
ft = zeros(NT,3);
ds = 2^(1-2*s)*gamma(1-s)/gamma(s);
if isfield(pde,'f') && ~ischar(pde.f)
    [lambda,weight] = quadpts(option.fquadorder);
    phif = lambda;                 % linear bases
    nQuadf = size(lambda,1);
    for p = 1:nQuadf
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for i = 1:3
            ft(:,i) = ft(:,i) + weight(p)*phif(p,i)*fp;
        end
    end
elseif isfield(pde,'f') && ischar(pde.f)
    if strcmp(pde.f,'intx') % f is a distribution 
        % compute (intx, d_x phi)
        [lambda,weight] = quadpts(option.fquadorder);
        nQuad = size(lambda,1);
        bt = zeros(NT,1);
        for p = 1:nQuad
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                + lambda(p,2)*node(elem(:,2),:) ...
                + lambda(p,3)*node(elem(:,3),:);
            bt = bt + weight(p)*pde.intx(pxy);      
        end
%         [Dphi,area] = gradbasis(node,elem);
        for i = 1:3
            ft(:,i) = Dphi(:,1,i).*bt; % (intx, d_xphi)
        end
    end
end
ft = ds*ft.*repmat(area,1,3);
b = b + accumarray(elem(:),ft(:),[Ndof,1]); 
% Neumann edfts are considered as open set. So the corner points should be
% set as Dirichlet boundary condition!
b(bdNode) = 0;
clear ft

%% Record assembling time
assembleTime = toc;  % stop of the timing
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
freeDof = find(bdidx==0);
if isempty(freeDof), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % Multigrid solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        tic;
        u(freeDof) = A(freeDof,freeDof)\b(freeDof);
        residual = norm(b - A*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        option.x0 = u;
        option.solver = 'VCYCLE';
        option.freeDof = freeDof;
        [u,info] = mgfracLapP1P1(A,b,elem,option);
    case 'amg'
        option.solver = 'CG';
        [u(freeNode),info] = amg(A(freeNode,freeNode),b(freeNode),option);                 
end

%% Compute error using boundary integral
if isfield(pde,'exactu')
    err = zeros(NT,1);
    [lambda,weight] = quadpts(7);
    phi = lambda;                 % linear bases
    nQuadf = size(lambda,1);
    for p = 1:nQuadf
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        uhp = u(elem(:,1))*phi(p,1) + u(elem(:,2))*phi(p,2) + u(elem(:,3))*phi(p,3);
        up = pde.exactu(pxy);
        err = err + weight(p)*fp.*(up - uhp);
%         err = err + weight(p)*fp.*(up - 2*uhp);
    end
    err = ds*sum(err.*area);
%     err = ds*sum(err.*area) + u'*(A*u);
    err = sqrt(err);
else
    err = 0;
end
info.errH1 = err;

%% Output information
eqn = struct('A',A,'b',b,'freeDof',freeDof,'y',y);
info.assembleTime = assembleTime;
info.aspectRatio = aspectRatio;