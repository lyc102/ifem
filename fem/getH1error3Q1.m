function err = getH1error3Q1(node,elem,Du,uh,K,option)
%% GETH1ERRORQ1 H1 norm of the approximation error on hex mesh.
%
%  err = getH1error3Q1(node,elem,@Du,uh) computes the H1 norm of the
%  error between the exact solution Du and finite element approximation
%  uh on a hex mesh described by node and elem. 
%
%  The input parameter Du is a function handle uh is a column array with
%  length N representing a piecewise linear function on the mesh.
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
% Modified by Long Chen. The original computation is not quite right. The
% diffusion coefficient cannot be computed separately.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Set up
if ~exist('option','var'), option = []; end
[N,  dim] = size(node);
[NT, nv] = size(elem);
if ~isfield('option','quadOrder')
   if exist('K','var') && ~isempty(K) && ~isnumeric(K)
        quadOrder = 3; 
   else
        quadOrder = 2; % number of Guass pts in one direction
   end
else
   quadOrder = option.quadorder;  
end

%% Mesh size and Volume
hx = node(elem(:,2),1) - node(elem(:,1),1);
hy = node(elem(:,4),2) - node(elem(:,1),2);
hz = node(elem(:,5),3) - node(elem(:,1),3);
volume = abs(hx.*hy.*hz);

%% Compute error
[GaussPt, weight] = quadptscube(quadOrder);
lambdax = GaussPt(:,1);
lambday = GaussPt(:,2);
lambdaz = GaussPt(:,3);

nQuad = size(GaussPt,1);
err = zeros(NT,1);
for p = 1:nQuad
    % Compute gradient of basis.
    if isfield(option,'isoparametric') && option.isoparametric
        [phi, Dphi, J] = hexbasis(node,elem,GaussPt(p,:));
    else
        Dphi = zeros(NT,3,8);
        phi2d = kron([lambday(p) 1-lambday(p)],[lambdax(p); 1-lambdax(p)]);
        phi = transpose(kron([lambdaz(p) 1-lambdaz(p)],phi2d([1 2 4 3])));
        Dx2d = kron([lambday(p) 1-lambday(p)],1./hx*[-1 1]);
        Dphi(:,1,:) = kron([lambdaz(p) 1-lambdaz(p)],Dx2d(:,[1 2 4 3]));
        Dy2d = kron(1./hy*[-1 1],[lambdax(p) 1-lambdax(p)]);
        Dphi(:,2,:) = kron([lambdaz(p) 1-lambdaz(p)],Dy2d(:,[1 2 4 3]));
        Dphi(:,3,:) = kron(1./hz*[-1 1],phi2d([1 2 4 3]));   
        J = volume;
    end    
%     [phi, Dphi, J] = hexbasis(node,elem, GaussPt(p,:));
    pxyz = zeros(NT, dim);
    for i = 1:dim
        xi = node(:,i);
        pxyz(:,i) = xi(elem)*phi;
    end
    Duh = zeros(NT,dim);
    for i = 1:dim
        Dphi_i(:,1:nv) = Dphi(:,i,:);
        Duh(:,i) = sum(uh(elem).*Dphi_i,2);
    end
    if exist('K','var') && ~isempty(K) && ~isnumeric(K) % K is a function
        err = err + weight(p)*K(pxyz).*sum((Du(pxyz)-Duh).^2,2);
    else
        err = err + weight(p)*sum((Du(pxyz)-Duh).^2,2);        
    end
end
if exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == NT
    err = K.*err;    % K is piecewise constant
end
err = J.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));