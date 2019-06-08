function [u,eqn,info] = fracLapP1P2(node,elem,pde,option)
%% FRACLAPP1P2 solves fractional Laplacian equation using P1-P2 element
%
% [u,eqn,info] = fracLapP1P1(node,elem,pde,option) solves the fractional
% Laplacian equation
% 
%  (-\Delta^^s u = f in \Omega with u = 0 on \partial \Omega.


%% Options
if ~exist('option','var'), option = []; end
if ~isfield(option,'plotflag'), option.plotflag = 0; end

%% Parameters
global s;
alpha = 1-2*s;
if s == 0.5
    Gamma = 1;
else
    Gamma = 5/(2*s)+0.01;  % a different grading factor
end

%% Mesh
N = size(node,1); NT = size(elem,1);
% My = round(sqrt(NT));
L = 1+log10(NT)/3; % to control the truncated error, L should increase as logN
My = round(L*2*NT^(1/4)); % number of elements in y direction
% My = round(pde.L*sqrt(sqrt(NT/2)));
y = gradmap(0,L,Gamma,My);
Ny = length(y); % number of vertices in y-direction  
NTy = Ny - 1;   % number of elements in y-direction
NTtotal = NT*NTy;
Nxdof = N;      % dof in x direction is number of vertices
Nydof = (Ny + NTy); % dof in y direction is number of vertices plus elements
Ndof = Nxdof*Nydof;

tic; % start the timing for the assembling
%% Stiffness matrix and mass matrix in the extended direction
% quantities in y direction
hy = diff(y);
a = zeros(NTy,5);
for i = 1:5
    a(:,i) = diff(y.^(alpha+i)/(alpha+i));
end
y1 = y(1:end-1);
y2 = y(2:end);
y12 = y1 + y2;
y1y2 = y1.*y2;
ym = y12/2;
% stiffness matrix in y direction
Ay = zeros(NTy,3,3);
Ay(:,1,1) = a(:,1)./hy.^2;
Ay(:,1,2) = -Ay(:,1,1);           
Ay(:,2,1) = -Ay(:,1,1);           
Ay(:,2,2) = Ay(:,1,1);         
Ay(:,1,3) = 8*(a(:,2) - ym.*a(:,1))./hy.^3;
Ay(:,3,1) = Ay(:,1,3);
Ay(:,2,3) = -Ay(:,1,3);
Ay(:,3,2) = Ay(:,2,3);
Ay(:,3,3) = 64*(a(:,3) - 2*ym.*a(:,2) + ym.^2.*a(:,1))./hy.^4;
% mass matrix in y direction
My = zeros(NTy,3,3);
My(:,1,1) = (a(:,3) - 2*y1.*a(:,2) + y1.^2.*a(:,1))./hy.^2;
My(:,1,2) = (-a(:,3) + (y1+y2).*a(:,2) - y1.*y2.*a(:,1))./hy.^2;
My(:,2,1) = My(:,1,2);
My(:,2,2) = (a(:,3) - 2*y2.*a(:,2) + y2.^2.*a(:,1))./hy.^2;
My(:,1,3) = 4*(a(:,4) - (y12+y2).*a(:,3) + (y2.*y12 + y1y2).*a(:,2) ...
             - y2.*y1y2.*a(:,1))./hy.^3;
My(:,3,1) = My(:,1,3);
My(:,2,3) = 4*(-a(:,4) + (y12+y1).*a(:,3) - (y1.*y12 + y1y2).*a(:,2) ...
             + y1.*y1y2.*a(:,1))./hy.^3;
My(:,3,2) = My(:,2,3);
My(:,3,3) = 16*(a(:,5) - 2*y12.*a(:,4) + (y12.^2 + 2*y1y2).*a(:,3)...
               -2*y12.*y1.*y2.*a(:,2) + y1y2.^2.*a(:,1))./hy.^4;

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

%% Assemble stiffness matrix
dofMap = repmat(1:Nxdof,Nydof,1) + repmat((0:Nydof-1)'*Nxdof,1,Nxdof);
ii = zeros(9*9*NTtotal,1); 
jj = zeros(9*9*NTtotal,1); 
sA = zeros(9*9*NTtotal,1);
index = 0;
for m = 1:3
    for n = 1:3
        for i = 1:3
            for j = 1:3
                if m < 3
                    ii(index+1:index+NTtotal) = dofMap(m:Ny+m-2,elem(:,i));
                else % m = 3: middle points of elements in y direction
                    ii(index+1:index+NTtotal) = dofMap(Ny+1:Nydof,elem(:,i));
                end
                if n < 3
                    jj(index+1:index+NTtotal) = dofMap(n:Ny+n-2,elem(:,j));
                else % n = 3: middle points of elements in y direction
                    jj(index+1:index+NTtotal) = dofMap(Ny+1:Nydof,elem(:,j));
                end                    
                sA(index+1:index+NTtotal) = kron(At(:,i,j),My(:,m,n)) + ...
                                            kron(Mt(:,i,j),Ay(:,m,n));
                index = index + NTtotal;
            end
        end
    end
end
A = sparse(ii,jj,sA,Ndof,Ndof);
clear ii jj sA At Mt Ay My

%% Assemble the right hand side
b = zeros(Ndof,1);
u = zeros(Ndof,1);

%% Set up boundary conditions
% find Dirichlet boundary nodes
bdNode = findboundary(elem);  % boundary vertices of mesh in x-direction
lateralbd = dofMap(:,bdNode);
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
    [lambdaf,weightf] = quadpts(option.fquadorder);
    phif = lambdaf;                 % linear bases
    nQuadf = size(lambdaf,1);
    for p = 1:nQuadf
        % quadrature points in the x-y coordinate
        pxy = lambdaf(p,1)*node(elem(:,1),:) ...
            + lambdaf(p,2)*node(elem(:,2),:) ...
            + lambdaf(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for i = 1:3
            ft(:,i) = ft(:,i) + weightf(p)*phif(p,i)*fp;
        end
    end
elseif isfield(pde,'f') && ischar(pde.f)
    if strcmp(pde.f,'intx') % f is a distribution 
        % compute (intx, d_x phi)
        [lambdaf,weightf] = quadpts(option.fquadorder);
        nQuad = size(lambdaf,1);
        bt = zeros(NT,1);
        for p = 1:nQuad
            pxy = lambdaf(p,1)*node(elem(:,1),:) ...
                + lambdaf(p,2)*node(elem(:,2),:) ...
                + lambdaf(p,3)*node(elem(:,3),:);
            bt = bt + weightf(p)*pde.intx(pxy);      
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
assembleTime = toc;
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
    else            % MGCG  solver for large size systems
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
        [u,info] = mgfracLapP1P2(A,b,elem,option); 
    case 'amg'
        option.solver = 'CG';
        [u(freeDof),info] = amg(A(freeDof,freeDof),b(freeDof),option);                 
end

%% Compute error using boundary integral
if isfield(pde,'exactu')
    err = zeros(NT,1);
    [lambda,weight] = quadpts(8);
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
    err = sqrt(abs(err));
else
    err = 0;
end
info.errH1 = err;

%% Output information
eqn = struct('A',A,'b',b,'freeDof',freeDof);
info.assembleTime = assembleTime;