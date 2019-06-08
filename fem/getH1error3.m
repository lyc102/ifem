function err = getH1error3(node,elem,Du,uh,K,quadOrder)
%% GETH1ERROR3 H1 norm of the approximation error in 3-D.
%
%  err = getH1error3(node,elem,@Du,uh) computes the H1 norm of the
%  error between the exact solution Du and finite element approximation
%  uh on a mesh described by node and elem. 
%
%  The input parameter Du is a function handle uh is a column array with
%  length N representing a piecewise linear function on the mesh or a
%  (NT,3) array representing the gradient of a linear function.
%
%  err = getH1error3(node,elem,@Du,uh,quadOrder) computes error
%  using the quadrature rule with order quadOrder (up to 5). The default
%  order is 2.
% 
%  Example: compute H1 error of piecewise linear interpolation
%
% [node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
% for k = 1:4
%     exactu = inline('sin(pi*p(:,1)).*sin(pi*p(:,2)).*p(:,3).^2','p');
%     Du = inline('[pi*cos(pi*p(:,1)).*sin(pi*p(:,2)).*p(:,3).^2, pi*sin(pi*p(:,1)).*cos(pi*p(:,2)).*p(:,3).^2, 2*sin(pi*p(:,1)).*sin(pi*p(:,2)).*p(:,3)]','p');
%     uI = exactu(node);
%     N(k) = size(node,1);
%     err(k) = getH1error3(node,elem,Du,uI);
%     [node,elem] = uniformbisect3(node,elem);
% end
% showrate(N,err);
%
% See also getH1error3, getL2error, getL2error3, quadpts.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

Nu = size(uh,1);    N = size(node,1);   NT = size(elem,1); 
% rough estimate using Euler formula N-NE+NF-NT = c, NF ~ 2NT, NE ~ N+NT-c
NE = N+NT;         NF = 2*NT;          NP2 = N + NE;  
if (Nu > 2*N+NT-5) && (Nu < 2*NT)    % Nu ~ N + NE: P2 element
    elem2dof = dof3P2(elem);
    NP2 = max(elem2dof(:));
end
if Nu > 2*NT % CR element or WG element
   elem2face = dof3face(elem);
   NF = max(elem2face(:));
end

%% Default quadrature orders for different elements
if ~exist('quadOrder','var')
    switch Nu
        case NT     % piecewise constant vector (uh is Duh)
            quadOrder = 2;
        case N      % piecewise linear function P1 element
            quadOrder = 2; 
        case NF     % piecewise linear function CR element
            quadOrder = 2; 
        case NP2    % piecewise quadratic function
            quadOrder = 3;
        case NF + NT % WG element
            quadOrder = 3;
    end
end

%% compute gradient of finite element function uh
if (size(uh,2) == 3) && (Nu == NT)      % uh is a piecewise constant vector
    Duh = uh;
    volume = abs(simplexvolume(node,elem));
elseif size(uh,2) == 1   % scalar function uh
    switch Nu
        case N      % piecewise linear function P1 element
            [Duh,volume] = gradu3(node,elem,uh);
        case NF     % piecewise linear function CR element
            [Duh,volume] = gradu3CR(node,elem,elem2face,uh); 
        case NF + NT % weak Galerkin element
            [Duh,volume] = gradu3WG(node,elem,elem2face,uh);             
        case NP2    % piecewise quadratic function
            [Dlambda,volume] = gradbasis3(node,elem);   
    end
end

%% compute H1 error element-wise using quadrature rule with order quadOrder
[lambda,weight] = quadpts3(quadOrder);
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
         + lambda(p,2)*node(elem(:,2),:) ...
         + lambda(p,3)*node(elem(:,3),:) ...
         + lambda(p,4)*node(elem(:,4),:);
    if Nu == NP2 % piecewise quadratic function
        Dphip1 = (4*lambda(p,1)-1).*Dlambda(:,:,1);
        Dphip2 = (4*lambda(p,2)-1).*Dlambda(:,:,2);
        Dphip3 = (4*lambda(p,3)-1).*Dlambda(:,:,3);
        Dphip4 = (4*lambda(p,4)-1).*Dlambda(:,:,4);
        Dphip5 = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
        Dphip6 = 4*(lambda(p,1)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,1));
        Dphip7 = 4*(lambda(p,1)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,1));
        Dphip8 = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
        Dphip9 = 4*(lambda(p,2)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,2));
        Dphip10 = 4*(lambda(p,3)*Dlambda(:,:,4)+lambda(p,4)*Dlambda(:,:,3));
        Duh = repmat(uh(elem2dof(:,1)),1,3).*Dphip1 + ...
              repmat(uh(elem2dof(:,2)),1,3).*Dphip2 + ...
              repmat(uh(elem2dof(:,3)),1,3).*Dphip3 + ...
              repmat(uh(elem2dof(:,4)),1,3).*Dphip4 + ...
              repmat(uh(elem2dof(:,5)),1,3).*Dphip5 + ...
              repmat(uh(elem2dof(:,6)),1,3).*Dphip6 + ...
              repmat(uh(elem2dof(:,7)),1,3).*Dphip7 + ...
              repmat(uh(elem2dof(:,8)),1,3).*Dphip8 + ...
              repmat(uh(elem2dof(:,9)),1,3).*Dphip9 + ...
              repmat(uh(elem2dof(:,10)),1,3).*Dphip10;        
    end     
    if exist('K','var') && ~isempty(K) && ~isreal(K) % K is a function
        err = err + weight(p)*K(pxyz).*sum((Du(pxyz)-Duh).^2,2);
    else
        err = err + weight(p)*sum((Du(pxyz)-Duh).^2,2);
    end
end
if exist('K','var') && ~isempty(K) && isreal(K) && size(K,1) == NT
    err = K.*err;    % K is piecewise constant
end
err = volume.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));