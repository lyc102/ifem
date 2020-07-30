function [err, errK] = getL2error3(node,elem,uexact,uh,quadOrder,varargin)
%% GETL2ERROR3 L2 norm of the approximation error in 3-D.
%
%  err = getL2error3(node,elem,@uexact,uh) computes the L2 norm of the error
%  between the exact solution uexact and a finite element approximation uh
%  on a mesh described by node and elem.
%  
%  [err, errK] = getL2error3(node,elem,@uexact,uh) also outputs the local error.
%
%  The input parameter uexact is a function handle and uh is a column array
%  which could be:
%    - P0 element: discontinuous and piecewise constant element
%    - P1 element: continuous linear element
%    - CR element: nonconforming linear element
%    - P2 element: continuous quadratic element
%    - P1+P0 element in Stokes equations
%
%  err = computeL2error3(node,elem,@uexact,uh,quadOrder) computes error
%  using the quadrature rule with order quadOrder (up to 5). The default
%  order is 3.
%   
%  Example: compute L2 error of nodal interpolation of linear element
%
%     exactu = inline('sin(pi*p(:,1)).*sin(pi*p(:,2)).*p(:,3).^2','p');
%     [node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
%     err = zeros(4,1); h = zeros(4,1);
%     for k = 1:4
%         [node,elem] = uniformrefine3(node,elem);
%         uI = exactu(node);
%         h(k) = size(elem,1)^(-1/3);
%         err(k) = getL2error3(node,elem,exactu,uI);
%     end
%     showrateh(h,err);
%
%  Vector Lagrange elements are also supported which is useful for Stokes
%  and elasticity equations.
%
%     pde = Maxwelldata3;
%     [node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
%     err = zeros(4,1); h = zeros(4,1);
%     for k = 1:4
%         [node,elem] = uniformrefine3(node,elem);
%         [elem2edge,edge] = dof3edge(elem);
%         uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
%         h(k) = size(elem,1)^(-1/3);
%         err(k) = getL2error3(node,elem,pde.exactu,uI);
%     end
%     showrateh(h,err);
%
% See also getL2error, getH1error, getH1error3, quadpts3.
%  
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

Nu = size(uh,1);    N = size(node,1);   NT = size(elem,1); 
% rough estimate using Euler formula 
% N-NE+NF-NT = c, NF = 4NT/2 + Nbf, NE = N + NT -c + Nbf
Nbf = round(NT^(2/3))-3; % estimate bd faces
NF = 2*NT + Nbf; 
NE = N + NT + Nbf - 1; 
NP2 = N + NE;  
if Nu > N % more than P1 element
   if Nu > 2*NT % element containing faces
       elem2face = dof3face(elem);
       NF = max(elem2face(:));
   end
   if  (Nu ~= NF) &&  (Nu ~= NF+NT) && (Nu ~= N+NT) && (Nu ~= NT)
       % try P2 element
        elem2dof = dof3P2(elem);
        NP2 = max(elem2dof(:));
   end
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
d = size(uh,2);  % uh could be a vector function
err = zeros(NT,d);
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
            uhp = uh(elem(:,1),:)*phi(p,1) + ...
                  uh(elem(:,2),:)*phi(p,2) + ...
                  uh(elem(:,3),:)*phi(p,3) + ...
                  uh(elem(:,4),:)*phi(p,4);
        case N + NT % P1+P0
            uhp = uh(elem(:,1),:)*phi(p,1) + ...
                  uh(elem(:,2),:)*phi(p,2) + ...
                  uh(elem(:,3),:)*phi(p,3) + ...
                  uh(elem(:,4),:)*phi(p,4);
            uhp = uhp + uh(N+1:end);        
        case NF  % CR nonconforming P1 element
            uhp = uh(elem2face(:,1),:)*phi(p,1) + ...
                  uh(elem2face(:,2),:)*phi(p,2) + ...
                  uh(elem2face(:,3),:)*phi(p,3) + ...
                  uh(elem2face(:,4),:)*phi(p,4);
        case NP2 % P2 piecewise quadratic function
            uhp = uh(elem2dof(:,1),:)*phi(p,1) + ...
                  uh(elem2dof(:,2),:)*phi(p,2) + ...
                  uh(elem2dof(:,3),:)*phi(p,3) + ...        
                  uh(elem2dof(:,4),:)*phi(p,4) + ...
                  uh(elem2dof(:,5),:)*phi(p,5) + ...
                  uh(elem2dof(:,6),:)*phi(p,6) + ...
                  uh(elem2dof(:,7),:)*phi(p,7) + ...                  
                  uh(elem2dof(:,8),:)*phi(p,8) + ...
                  uh(elem2dof(:,9),:)*phi(p,9) + ...                  
                  uh(elem2dof(:,10),:)*phi(p,10);
        case NF+NT  % weak Galerkin element
%             uhp = uh(1:NT); % only count the interior part
            uhp = uh(elem2face(:,1),:)*phi(p,1) + ...
                  uh(elem2face(:,2),:)*phi(p,2) + ...
                  uh(elem2face(:,3),:)*phi(p,3) + ...
                  uh(elem2face(:,4),:)*phi(p,4);              
    end
    % quadrature points in the x-y coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
         + lambda(p,2)*node(elem(:,2),:) ...
         + lambda(p,3)*node(elem(:,3),:) ...
         + lambda(p,4)*node(elem(:,4),:);
    % when size(uhp,2) == 1, summing on 2-nd component has no impact
    err = err + weight(p)*sum((uexact(pxyz) - uhp).^2,2);
end

%% Modification
% volume of tetrahedrons
d12 = node(elem(:,2),:)-node(elem(:,1),:);
d13 = node(elem(:,3),:)-node(elem(:,1),:);
d14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = abs(dot(mycross(d12,d13,2),d14,2))/6;
err = volume.*sum(err,2);
err(isnan(err)) = 0; % singular point is excluded
if nargout > 1; errK = sqrt(err); end
err = sqrt(sum(err));
