function err = getL2error(node,elem,uexact,uh,quadOrder,varargin)
%% GETL2ERROR L2 norm of the approximation error.
%
%  err = getL2error(node,elem,@uexact,uh) computes the L2 norm of the error
%  between the exact solution uexact and a finite element approximation uh
%  on a mesh described by node and elem.
%
%  The input parameter uexact is a function handle and uh is a column array
%  which could be:
%    - P0 element i.e. discontinuous and piecewise constant
%    - P1 element i.e. continuous and piecewise linear
%    - CR element i.e. piecewise linear and continuous at mid pts of edges
%    - P2 element i.e. continuous and piecewise quadratic
%    - P1+P0 element which could happen in the fluid application
%
%  err = getL2error(node,elem,@uexact,uh,quadOrder) computes error
%  using the quadrature rule with order quadOrder (up to 5). The default
%  order is 3.
%   
%  Example: compute L2 error of piecewise linear interpolation
%
%     [node,elem] = squaremesh([0,1,0,1],0.25);
%     for k = 1:4
%         exactu = inline('sin(pi*pxy(:,1)).*sin(pi*pxy(:,2))','pxy');
%         uI = exactu(node);
%         N(k) = size(node,1);
%         err(k) = getL2error(node,elem,exactu,uI);
%         [node,elem] = uniformrefine(node,elem);
%     end
%     showrate(N,err);
%
% The cubic element is added by Jie Zhou.
%
% See also getH1error, getL2error3, getH1error3, quadpts.
%  
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Number of vertices, elements, edges, and degrees of freedom
Nu = length(uh);    N = size(node,1);   NT = size(elem,1);   
% Euler formula N-NE+NT = c % rough estimateus using Euler formula
NE = N + NT - 1;    NP2 = N + NE;   NP3 = N + 2*NE + NT;  
if Nu > N+NT-5
    elem2dof = dofP2(elem);
    NP2 = max(elem2dof(:));
    NE = NP2 - N;
    NP3 = N+2*NE+NT;   
end

%% Default quadrature orders for different elements
if ~exist('quadOrder','var') || isempty(quadOrder)
    switch Nu
        case NT     % piecewise constant function P0
            quadOrder = 2;
        case N      % piecewise linear function P1 element
            quadOrder = 3; 
        case NE     % piecewise linear function CR element
            quadOrder = 3; 
        case N+NT   % piecewise linear function + constant function
            quadOrder = 3;        
        case NP2    % piecewise quadratic function
            quadOrder = 4;
        case NE+NT  % weak Galerkin element
            quadOrder = 3;
        case NP3    % P3 element
            quadOrder = 5;            
    end
end

%% compute L2 error element-wise using quadrature rule with order quadOrder
err = zeros(NT,1);
[lambda,weight] = quadpts(quadOrder);
% basis function at quadrature points
switch Nu
    case N    % P1 piecewise linear function
        phi = lambda; % linear bases
    case N+NT % P1+P0
        phi = lambda; % linear bases
    case NE  % CR nonconforming P1 element
        phi = 1-2*lambda;
        elem2edge = elem2dof(:,4:6) - N;
    case NP2 % P2 piecewise quadratic elements
        phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
        phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
        phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
        phi(:,4) = 4*lambda(:,2).*lambda(:,3);
        phi(:,5) = 4*lambda(:,1).*lambda(:,3);
        phi(:,6) = 4*lambda(:,2).*lambda(:,1);
    case NE+NT  % weak Galerkin element
%             uhp = uh(1:NT); % only count the interior part
        phi = 1-2*lambda;
        elem2edge = elem2dof(:,4:6) - N + NT;        
    case 2*NE+NT+N % P3 piecewise cubic elements          
        phi(:,1) = 0.5*(3*lambda(:,1)-1).*(3*lambda(:,1)-2).*lambda(:,1);           
        phi(:,2) = 0.5*(3*lambda(:,2)-1).*(3*lambda(:,2)-2).*lambda(:,2); 
        phi(:,3) = 0.5*(3*lambda(:,3)-1).*(3*lambda(:,3)-2).*lambda(:,3);
        phi(:,4) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,2)-1); 
        phi(:,5) = 9/2*lambda(:,3).*lambda(:,2).*(3*lambda(:,3)-1); 
        phi(:,6) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,3)-1);      
        phi(:,7) = 9/2*lambda(:,1).*lambda(:,3).*(3*lambda(:,1)-1);  
        phi(:,8) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,1)-1);
        phi(:,9) = 9/2*lambda(:,1).*lambda(:,2).*(3*lambda(:,2)-1);        
        phi(:,10) = 27*lambda(:,1).*lambda(:,2).*lambda(:,3);    
        elem2dof = dofP3(elem);           
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
                  uh(elem(:,3))*phi(p,3);
        case N+NT % P1+P0
            uhp = uh(elem(:,1))*phi(p,1) + ...
                  uh(elem(:,2))*phi(p,2) + ...
                  uh(elem(:,3))*phi(p,3);
            uhp = uhp + uh(N+1:end);
        case NE  % CR nonconforming P1 element
            uhp = uh(elem2edge(:,1))*phi(p,1) + ...
                  uh(elem2edge(:,2))*phi(p,2) + ...
                  uh(elem2edge(:,3))*phi(p,3);
        case NP2 % P2 piecewise quadratic function
            uhp = uh(elem2dof(:,1)).*phi(p,1) + ...
                  uh(elem2dof(:,2)).*phi(p,2) + ...
                  uh(elem2dof(:,3)).*phi(p,3) + ...        
                  uh(elem2dof(:,4)).*phi(p,4) + ...
                  uh(elem2dof(:,5)).*phi(p,5) + ...
                  uh(elem2dof(:,6)).*phi(p,6);
        case NP3
            uhp = uh(elem2dof(:,1)).*phi(p,1) + ...
                  uh(elem2dof(:,2)).*phi(p,2) + ...
                  uh(elem2dof(:,3)).*phi(p,3) + ...
                  uh(elem2dof(:,4)).*phi(p,4) + ...
                  uh(elem2dof(:,5)).*phi(p,5) + ...
                  uh(elem2dof(:,6)).*phi(p,6) + ...
                  uh(elem2dof(:,7)).*phi(p,7) + ...
                  uh(elem2dof(:,8)).*phi(p,8) + ...
                  uh(elem2dof(:,9)).*phi(p,9) + ...
                  uh(elem2dof(:,10)).*phi(p,10);
        case NE+NT  % weak Galerkin element
%             uhp = uh(1:NT); % only count the interior part
            uhp = uh(elem2edge(:,1))*phi(p,1) + ...
                  uh(elem2edge(:,2))*phi(p,2) + ...
                  uh(elem2edge(:,3))*phi(p,3);
    end
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    err = err + weight(p)*(uexact(pxy) - uhp).^2;
end
%% Modification
% area of triangles
ve2 = node(elem(:,1),:)-node(elem(:,3),:);
ve3 = node(elem(:,2),:)-node(elem(:,1),:);
area = 0.5*abs(-ve3(:,1).*ve2(:,2)+ve3(:,2).*ve2(:,1));
err = area.*err;
err(isnan(err)) = 0; % singular values, i.e. uexact(p) = infty, are excluded
err = sqrt(sum(err));