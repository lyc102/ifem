function [eta,estar,oscf] = estimatefracLap(node,elem,y,u,pde,option)

if ~exist('option','var'), option = []; end

%% Parameters
global s
alpha = 1-2*s;

%% Mesh and data structure
[elem2dof,edge,bdDof] = dofP2(elem);
N = size(node,1);  NT = size(elem,1); NE = size(edge,1);
Ny = length(y); % number of vertices in y-direction  
NTy = Ny - 1;   % number of elements in y-direction
NTtotal = NT*NTy;
Nxdof = (N + NE + NT);
Nydof = (Ny + NTy);
Ndof = Nxdof*Nydof;
NTdof = N + NE + (1:NT)';
elem2dof(:,7) = NTdof; % add elementwise bubble

tic;
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
[Dlambda,area] = gradbasis(node,elem);
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
if ~isempty(pde.d) % build the coefficients into the area
    areaK = K.*area;
else
    areaK = area;
end
At = zeros(NT,7,7);
Mt = zeros(NT,7,7);
% stiffness matrix
[lambda, w] = quadpts(4);
nQuad = size(lambda,1);
for p = 1:nQuad
    % Dphi at quadrature points
    Dphip(:,:,7) = 27*(lambda(p,1)*lambda(p,2)*Dlambda(:,:,3) + ...
                       lambda(p,1)*lambda(p,3)*Dlambda(:,:,2) + ...
                       lambda(p,3)*lambda(p,2)*Dlambda(:,:,1));           
    Dphip(:,:,1) = (4*lambda(p,1)-1).*Dlambda(:,:,1);            
    Dphip(:,:,2) = (4*lambda(p,2)-1).*Dlambda(:,:,2);            
    Dphip(:,:,3) = (4*lambda(p,3)-1).*Dlambda(:,:,3);            
    Dphip(:,:,4) = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
    Dphip(:,:,5) = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
    Dphip(:,:,6) = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
    for i = 1:7
        for j = 1:7
            At(:,i,j) = At(:,i,j) + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*areaK;           
        end
    end
end
% mass matrix
[lambda, w] = quadpts(6);
nQuad = size(lambda,1);
phi(:,7) = 27*lambda(:,1).*lambda(:,2).*lambda(:,3); 
phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
phi(:,4) = 4*lambda(:,2).*lambda(:,3);
phi(:,5) = 4*lambda(:,3).*lambda(:,1);
phi(:,6) = 4*lambda(:,1).*lambda(:,2);
for i = 1:7
    for j = 1:7
        for p = 1:nQuad
            Mt(:,i,j) = Mt(:,i,j) + w(p)*phi(p,i).*phi(p,j);
        end
        Mt(:,i,j) = Mt(:,i,j).*areaK;
    end
end
clear phi Dphip lambda Dlambda

%% Assemble stiffness matrix
dofMap = repmat(1:Nxdof,Nydof,1) + repmat((0:Nydof-1)'*Nxdof,1,Nxdof);
ii = zeros(21*21*NTtotal,1); 
jj = zeros(21*21*NTtotal,1); 
sA = zeros(21*21*NTtotal,1);
index = 0;
for m = 1:3
    for n = 1:3
        for i = 1:7
            for j = 1:7
                if m < 3
                    ii(index+1:index+NTtotal) = dofMap(m:Ny+m-2,elem2dof(:,i));
                else % m = 3: middle points of elements in y direction
                    ii(index+1:index+NTtotal) = dofMap(Ny+1:Nydof,elem2dof(:,i));
                end
                if n < 3
                    jj(index+1:index+NTtotal) = dofMap(n:Ny+n-2,elem2dof(:,j));
                else % n = 3: middle points of elements in y direction
                    jj(index+1:index+NTtotal) = dofMap(Ny+1:Nydof,elem2dof(:,j));
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

%% Compute boundary integral over the original domain
if ~isfield(option,'fquadorder')
    option.fquadorder = 7;   
end
ft = zeros(NT,7);
ds = 2^(1-2*s)*gamma(1-s)/gamma(s);
if isfield(pde,'f') && ~ischar(pde.f)
    [lambda,weight] = quadpts(option.fquadorder);
    phif(:,7) = 27*lambda(:,1).*lambda(:,2).*lambda(:,3);
    phif(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phif(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phif(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phif(:,4) = 4*lambda(:,2).*lambda(:,3);
    phif(:,5) = 4*lambda(:,3).*lambda(:,1);
    phif(:,6) = 4*lambda(:,1).*lambda(:,2);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for i = 1:7
            ft(:,i) = ft(:,i) + weight(p)*phif(p,i)*fp;
        end
    end
elseif isfield(pde,'f') && ischar(pde.f)
    if strcmp(pde.f,'intx') % f is a distribution 
        % compute (intx, d_x phi)
        [lambda,weight] = quadpts(option.fquadorder);
        nQuad = size(lambda,1);
        [Dlambda,area] = gradbasis(node,elem);
        for p = 1:nQuad
            % Dphi at quadrature points
            Dxphip(:,7) = 27*(lambda(p,1)*lambda(p,2)*Dlambda(:,1,3) + ...
                              lambda(p,1)*lambda(p,3)*Dlambda(:,1,2) + ...
                              lambda(p,3)*lambda(p,2)*Dlambda(:,1,1));                       
            Dxphip(:,6) = 4*(lambda(p,1)*Dlambda(:,1,2)+lambda(p,2)*Dlambda(:,1,1));            
            Dxphip(:,1) = (4*lambda(p,1)-1).*Dlambda(:,1,1);            
            Dxphip(:,2) = (4*lambda(p,2)-1).*Dlambda(:,1,2);            
            Dxphip(:,3) = (4*lambda(p,3)-1).*Dlambda(:,1,3);            
            Dxphip(:,4) = 4*(lambda(p,2)*Dlambda(:,1,3)+lambda(p,3)*Dlambda(:,1,2));
            Dxphip(:,5) = 4*(lambda(p,3)*Dlambda(:,1,1)+lambda(p,1)*Dlambda(:,1,3));
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                + lambda(p,2)*node(elem(:,2),:) ...
                + lambda(p,3)*node(elem(:,3),:);            
            for i = 1:7
                ft(:,i) = ft(:,i) + weight(p)*Dxphip(:,i).*pde.intx(pxy);
            end
        end
    end
end
ft = ds*ft.*repmat(area,1,7);
b = b + accumarray(elem2dof(:),ft(:),[Ndof,1]); 
% Neumann edges are considered as open set. So the corner points should be
% set as Dirichlet boundary condition!
% b(lateralbd) = 0;
clear ft

%% Prolongate P1 solution to P2 space
u = reshape(u,N,Ny);
P1toP2 = sparse([(1:N)'; N+(1:NE)'; N+(1:NE)'], ...
                [(1:N)'; double(edge(:))],...
                [ones(N,1); 0.5*ones(2*NE,1)]',Nxdof,N);
u2 = P1toP2*u;
% the bubble function is added as a hierarchical bases. so no prolongation.
u2 = u2(:);
u2(Ndof) = 0; % prolongation in y-direction. zero extension since HB basis

%% Record assembling time
assembleTime = toc;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Local problems on stars
r = b - A*u2;  % form residual
isBdEdge = false(NE,1);  
isBdEdge(bdDof(bdDof>N)-N) = true; % find boundary edges
isBdNode = false(N,1);
isBdNode(bdDof(bdDof<=N)) = true;  % find boundary nodes
intEdgeIdx = find(~isBdEdge);
e2v = sparse([intEdgeIdx,intEdgeIdx], double(edge(intEdgeIdx,:)), 1, NE, N);
t2v = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), 1, NT, N);
estar = zeros(N,1);
for i = 1:N
    starEdges = find(e2v(:,i)); % include interior edges only
    starElem = find(t2v(:,i));
    if isBdNode(i) % for bd nodes, it is not included in the local dof
        starDof = dofMap([1:Ny-1,Ny+1:Nydof],[starEdges+N; starElem+N+NE]);
    else
        starDof = dofMap([1:Ny-1,Ny+1:Nydof],[i; starEdges+N; starElem+N+NE]);
    end
    localA = A(starDof(:),starDof(:));
    etastar = localA\r(starDof(:));
    estar(i) = etastar'*localA*etastar;
end

%% Data oscillation
oscf = zeros(NT,1);
if isfield(pde,'f') && ~ischar(pde.f) % exclude the distribution case
    % compute average
    ft = zeros(NT,1);
    if ~isfield(option,'dataquadorder')
        option.dataquadorder = 7;   
    end
    [lambda,weight] = quadpts(option.dataquadorder);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        ppxy = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ...
             + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(ppxy);
        ft = ft + weight(p)*fp;
    end
    % compute h^{2s}||f - ft||^2
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        ppxy = lambda(p,1)*node(elem(:,1),:) ...
             + lambda(p,2)*node(elem(:,2),:) ...
             + lambda(p,3)*node(elem(:,3),:)   ;
        fp = pde.f(ppxy);
        oscf = oscf + weight(p)*(fp - ft).^2;
    end
    oscf = ds*oscf.*area;
    h = sqrt(area);
    oscf = oscf.*h.^(2*s);
end

%% Change nodal wise to element wise indicator
valence = accumarray(elem(:),ones(3*NT,1),[N 1]);
estarT = estar./valence;
eta = estarT(elem(:,1)) + estarT(elem(:,2)) + estarT(elem(:,3));
% figure(3); hold off; plot(eta); hold on; plot(oscf,'-r');
% legend('Estar','osc')
eta = sqrt(eta + oscf);
estar = sqrt(estar);