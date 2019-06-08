function err = getL2error3(node,elem,uexact,uh,quadOrder,varargin)
%% GETL2ERROR3 L2 norm of the approximation error in 3-D.
%
%  err = getL2error3(node,elem,@uexact,uh) computes the L2 norm of the error
%  between the exact solution uexact and a finite element approximation uh
%  on a mesh described by node and elem.
%
%  The input parameter uexact is a function handle and uh is a column array
%  which could be:
%    - P0 element i.e. discontinuous and piecewise constant
%    - P1 element i.e. continuous and piecewise linear
%    - P2 element i.e. continuous and piecewise quadratic
%    - P1+P0 element which could happen in the fluid application
%
%  err = computeL2error3(node,elem,@uexact,uh,quadOrder) computes error
%  using the quadrature rule with order quadOrder (up to 5). The default
%  order is 3.
%   
%  Example: compute L2 error of piecewise linear interpolation
%
%     [node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
%     for k = 1:4
%         exactu = inline('sin(pi*p(:,1)).*sin(pi*p(:,2)).*p(:,3).^2','p');
%         uI = exactu(node);
%         N(k) = size(node,1);
%         err(k) = getL2error3(node,elem,exactu,uI);
%         [node,elem] = uniformrefine3(node,elem);
%     end
%     showrate(N,err);
%
% See also getL2error, getH1error, getH1error3, quadpts3.
%  
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

Nu = length(uh);    N = size(node,1);   NT = size(elem,1); 
% Euler formula N-NE+NF-NT = c, NF ~ 2NT, NE ~ N+NT-c
NE = N+NT; NF = 2*NT; NP2 = N + NE;  % rough estimate
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
        case NT     % piecewise constant function
            quadOrder = 2;
        case N      % piecewise linear function
            quadOrder = 3; 
        case NF     % piecewise linear function CR element
            quadOrder = 3;             
        case N+NT   % piecewise linear function + constant function
            quadOrder = 3;        
        case NP2    % piecewise quadratic function
            quadOrder = 4;
        case NF+NT  % WG element
            quadOrder = 3;      
    end
end
%% compute L2 error element-wise using quadrature rule with order quadOrder
err = zeros(NT,1);
[lambda,weight] = quadpts3(quadOrder);
% basis function at quadrature points
switch Nu
    case N    % P1 piecewise linear function
        phi = lambda; % linear bases
    case N+NT % P1+P0
        phi = lambda; % linear bases
    case NF  % CR nonconforming P1 element
        phi = 1-3*lambda;
    case NP2 % P2 piecewise quadratic function
        phi(:,10)= 4*lambda(:,3).*lambda(:,4);
        phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
        phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
        phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
        phi(:,4) = lambda(:,4).*(2*lambda(:,4)-1);
        phi(:,5) = 4*lambda(:,1).*lambda(:,2);
        phi(:,6) = 4*lambda(:,1).*lambda(:,3);
        phi(:,7) = 4*lambda(:,1).*lambda(:,4);
        phi(:,8) = 4*lambda(:,2).*lambda(:,3);
        phi(:,9) = 4*lambda(:,2).*lambda(:,4);
    case NF+NT  % weak Galerkin element
%             uhp = uh(1:NT); % only count the interior part
        phi = 1-3*lambda;
        elem2face = NT + elem2face;            
end
nQuad = size(lambda,1);
for p = 1:nQuad
    % evaluate uh at quadrature point
    switch Nu
        case NT   % P0 piecewise constant function
            uhp = uh;
        case N    % P1 piecewise linear function
            uhp = uh(elem(:,1))*phi(p,1) + ...
                  uh(elem(:,2))*phi(p,2) + ...
                  uh(elem(:,3))*phi(p,3) + ...
                  uh(elem(:,4))*phi(p,4);
        case N + NT % P1+P0
            uhp = uh(elem(:,1))*phi(p,1) + ...
                  uh(elem(:,2))*phi(p,2) + ...
                  uh(elem(:,3))*phi(p,3) + ...
                  uh(elem(:,4))*phi(p,4);
            uhp = uhp + uh(N+1:end);        
        case NF  % CR nonconforming P1 element
            uhp = uh(elem2face(:,1))*phi(p,1) + ...
                  uh(elem2face(:,2))*phi(p,2) + ...
                  uh(elem2face(:,3))*phi(p,3) + ...
                  uh(elem2face(:,4))*phi(p,4);
        case NP2 % P2 piecewise quadratic function
            uhp = uh(elem2dof(:,1))*phi(p,1) + ...
                  uh(elem2dof(:,2))*phi(p,2) + ...
                  uh(elem2dof(:,3))*phi(p,3) + ...        
                  uh(elem2dof(:,4))*phi(p,4) + ...
                  uh(elem2dof(:,5))*phi(p,5) + ...
                  uh(elem2dof(:,6))*phi(p,6) + ...
                  uh(elem2dof(:,7))*phi(p,7) + ...                  
                  uh(elem2dof(:,8))*phi(p,8) + ...
                  uh(elem2dof(:,9))*phi(p,9) + ...                  
                  uh(elem2dof(:,10))*phi(p,10);
        case NF+NT  % weak Galerkin element
%             uhp = uh(1:NT); % only count the interior part
            uhp = uh(elem2face(:,1))*phi(p,1) + ...
                  uh(elem2face(:,2))*phi(p,2) + ...
                  uh(elem2face(:,3))*phi(p,3) + ...
                  uh(elem2face(:,4))*phi(p,4);              
    end
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    err = err + weight(p)*(uexact(pxy) - uhp).^2;
end
%% Modification
% volume of tetrahedrons
d12 = node(elem(:,2),:)-node(elem(:,1),:);
d13 = node(elem(:,3),:)-node(elem(:,1),:);
d14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = abs(dot(mycross(d12,d13,2),d14,2)/6);
err(isnan(err)) = 0; % singular point is excluded
err = sqrt(sum(volume.*err));