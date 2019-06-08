function err = getH1error(node,elem,Du,uh,K,quadOrder)
%% GETH1ERROR H1 norm of the approximation error.
%
%  err = getH1error(node,elem,@Du,uh,K) computes the H1 norm of the
%  error between the exact solution Du and finite element approximation
%  uh on a mesh described by node and elem. 
%
%  The input parameter Du is a function handle and uh is either a column
%  array with length N representing a piecewise linear function on the mesh
%  or a (NT,2) array representing the gradient of a linear function. The
%  diffusion coefficient K is either a NT by 1 array or a function handle.
%
%  err = getH1error(node,elem,@Du,uh,d,quadOrder) computes error
%  using the quadrature rule with order quadOrder (up to 9). The default
%  order is 3.
% 
%  Example: compute H1 error of piecewise linear interpolation
%
%     [node,elem] = squaremesh([0,1,0,1],0.25);
%     for k = 1:4
%         exactu = inline('sin(pi*pxy(:,1)).*sin(pi*pxy(:,2))','pxy');
%         Du = inline('[pi*cos(pi*pxy(:,1)).*sin(pi*pxy(:,2)) pi*sin(pi*pxy(:,1)).*cos(pi*pxy(:,2))]','pxy');
%         uI = exactu(node);
%         N(k) = size(node,1);
%         err(k) = getH1error(node,elem,Du,uI);
%         [node,elem] = uniformrefine(node,elem);
%     end
%     showrate(N,err);
%
% See also getH1error3, getL2error, getL2error3, quadpts.
%
% The quadratic element is added by Ming Wang. The cubic element is added by Jie Zhou.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

Nu = size(uh,1);    N = size(node,1);   NT = size(elem,1);
% Euler formula N-NE+NT = c % rough estimateus using Euler formula
NE = N + NT;    NP2 = N + NE;   NP3 = N + 2*NE + NT;    
if Nu > N+NT-5   % Euler formula N-NE+NT = c
    elem2dof = dofP2(elem);
    NP2 = max(elem2dof(:));
    NE = NP2 - N;
    NP3 = N+2*NE+NT;    
end

%% Default quadrature orders for different elements
if ~exist('quadOrder','var')
    switch Nu
        case NT     % piecewise constant vector (uh is Duh)
            quadOrder = 3;
        case N      % piecewise linear function P1 element
            quadOrder = 3; 
        case NE     % piecewise linear function CR element
            quadOrder = 3; 
        case NP2    % piecewise quadratic function
            quadOrder = 4;
        case NE + NT % WG element
            quadOrder = 3;       
        case NP3    % P3 element
            quadOrder = 5;               
    end
end

%% compute gradient of finite element function uh
if (size(uh,2) == 2) && (Nu == NT)      % uh is a piecewise constant vector
    Duh = uh;
    area = abs(simplexvolume(node,elem));
elseif size(uh,2) == 1   % scalar function uh
    switch Nu
        case N      % piecewise linear function P1 element
            [Duh,area] = gradu(node,elem,uh);
        case NE     % piecewise linear function CR element
            elem2edge = elem2dof(:,4:6) - N;            
            [Duh,area] = graduCR(node,elem,elem2edge,uh); 
        case NE + NT % weak Galerkin element
            elem2edge = elem2dof(:,4:6) - N;            
            [Duh,area] = graduWG(node,elem,elem2edge,uh);             
        case NP2    % piecewise quadratic function
            [Dlambda,area] = gradbasis(node,elem);
        case NP3
            [Dlambda,area] = gradbasis(node,elem);   
            elem2dof = dofP3(elem);
    end
end

%% compute H1 error element-wise using quadrature rule with order quadOrder
[lambda,weight] = quadpts(quadOrder);
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    if Nu == NP2 % piecewise quadratic function
        Dphip1 = (4*lambda(p,1)-1).*Dlambda(:,:,1);
        Dphip2 = (4*lambda(p,2)-1).*Dlambda(:,:,2);
        Dphip3 = (4*lambda(p,3)-1).*Dlambda(:,:,3);
        Dphip4 = 4*(lambda(p,2)*Dlambda(:,:,3)+lambda(p,3)*Dlambda(:,:,2));
        Dphip5 = 4*(lambda(p,3)*Dlambda(:,:,1)+lambda(p,1)*Dlambda(:,:,3));
        Dphip6 = 4*(lambda(p,1)*Dlambda(:,:,2)+lambda(p,2)*Dlambda(:,:,1));
        Duh = repmat(uh(elem2dof(:,1)),1,2).*Dphip1 + ...
              repmat(uh(elem2dof(:,2)),1,2).*Dphip2 + ...
              repmat(uh(elem2dof(:,3)),1,2).*Dphip3 + ...
              repmat(uh(elem2dof(:,4)),1,2).*Dphip4 + ...
              repmat(uh(elem2dof(:,5)),1,2).*Dphip5 + ...
              repmat(uh(elem2dof(:,6)),1,2).*Dphip6;
    end
    if Nu == NP3 % piecewise cubic function
        Dphip1 = (27/2*lambda(p,1)*lambda(p,1)-9*lambda(p,1)+1).*Dlambda(:,:,1);           
        Dphip2 = (27/2*lambda(p,2)*lambda(p,2)-9*lambda(p,2)+1).*Dlambda(:,:,2); 
        Dphip3 = (27/2*lambda(p,3)*lambda(p,3)-9*lambda(p,3)+1).*Dlambda(:,:,3);
        Dphip4 = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda(:,:,3)+...
                lambda(p,3)*(6*lambda(p,2)-1).*Dlambda(:,:,2));  
        Dphip5 = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda(:,:,2)+...
                 lambda(p,2)*(6*lambda(p,3)-1).*Dlambda(:,:,3));             
        Dphip6 = 9/2*((3*lambda(p,3)*lambda(p,3)-lambda(p,3)).*Dlambda(:,:,1)+...
                 lambda(p,1)*(6*lambda(p,3)-1).*Dlambda(:,:,3)); 
        Dphip7 = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda(:,:,3)+...
                 lambda(p,3)*(6*lambda(p,1)-1).*Dlambda(:,:,1)); 
        Dphip8 = 9/2*((3*lambda(p,1)*lambda(p,1)-lambda(p,1)).*Dlambda(:,:,2)+...
                 lambda(p,2)*(6*lambda(p,1)-1).*Dlambda(:,:,1)); 
        Dphip9 = 9/2*((3*lambda(p,2)*lambda(p,2)-lambda(p,2)).*Dlambda(:,:,1)+...
                 lambda(p,1)*(6*lambda(p,2)-1).*Dlambda(:,:,2));  
        Dphip10= 27*(lambda(p,1)*lambda(p,2)*Dlambda(:,:,3)+lambda(p,1)*lambda(p,3)*Dlambda(:,:,2)+...
                 lambda(p,3)*lambda(p,2)*Dlambda(:,:,1));
        Duh = repmat(uh(elem2dof(:,1)),1,2).*Dphip1 + ...
              repmat(uh(elem2dof(:,2)),1,2).*Dphip2 + ...
              repmat(uh(elem2dof(:,3)),1,2).*Dphip3 + ...
              repmat(uh(elem2dof(:,4)),1,2).*Dphip4 + ...
              repmat(uh(elem2dof(:,5)),1,2).*Dphip5 + ...
              repmat(uh(elem2dof(:,6)),1,2).*Dphip6 + ...
              repmat(uh(elem2dof(:,7)),1,2).*Dphip7 + ...
              repmat(uh(elem2dof(:,8)),1,2).*Dphip8 + ...
              repmat(uh(elem2dof(:,9)),1,2).*Dphip9 + ...    
              repmat(uh(elem2dof(:,10)),1,2).*Dphip10;
    end 
    if exist('K','var') && ~isempty(K) && ~isnumeric(K) % K is a function
        err = err + weight(p)*K(pxy).*sum((Du(pxy)-Duh).^2,2);
    else
        err = err + weight(p)*sum((Du(pxy)-Duh).^2,2);        
    end
end
if exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == NT
    err = K.*err;    % K is piecewise constant
end
err = area.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));