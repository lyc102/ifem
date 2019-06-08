function eta = estimateresidual(node,elem,u,pde,bdFlag)
%% ESTIMATERESIDUAL residual type error estimator.
%
% eta = estimateresidual(node,elem,u,pde) computes the residual type a
% posterior error estimator for the Poisson equation in the form
%
% $$\eta _T^2 = \sum _{E\in T}[\nabla u\cdot v_E^{\bot}]^2 + area(T)^2\sum _i w_if_i^2,$$
%
% which is an approximation of 
%
% $$\eta _T^2 =\sum _{E\in T}\int _E h_E[\nabla u\cdot n_E]^2 ds + \int _T h^2f^2 dx.$$
%
% Here [\nabla u\cdot n_E] denotes the jump of the flux accross the
% interior edge E and n_E is the normal vector of E. In 2-D, n_E is the
% right 90 degree rotation of the edge vector of E. The integral $\int _E
% h_E[\nabla u\cdot n_E]^2 ds$ is simplified to $[\nabla u\cdot
% v_E^{\bot}]^2$.
%
% On Neumann boundary edges, the jump as is modfied to g - du/dn and the
% contribution from the Neumann edge is area(T)*|g - du/dn|^2. In 3-D,
% h_E^2 is approximated by area(T).
%
% For general diffusion equation div(K grad u) = f, the flux is K grad u.
% The jump will be [K\nabla u\cdot n_E].
%
% PDE is a structure array that records the information of Neumann
% boundary condition and diffusion coefficent.
%
%  pde.f     : right hand side
%  pde.g_N   : Neumann boundary condition data
%  pde.d     : diffusion coeffcients
%
% Example
%
%   Kellogg
%
% See also estimateresidual3, estimaterecovery, estimaterecovery3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Compute the flux from the finite element function u
[Du,area] = gradu(node,elem,u);
if exist('pde','var') && isfield(pde,'d') && ~isempty(pde.d)
    if isreal(pde.d)
        K = pde.d;                  % d is an array
    else                            % d is a function
        center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
        K = pde.d(center);              
    end
    Du = [K.*Du(:,1),K.*Du(:,2)];   % flux
end

%% Jump of normal flux
% data structure
T = auxstructure(elem);
neighbor = T.neighbor; 
clear T
% edge vector
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
% scaled normal vector: a right 90 degree rotation of edge vector
ne(:,1,1)= ve(:,2,1); ne(:,2,1)= -ve(:,1,1);
ne(:,1,2)= ve(:,2,2); ne(:,2,2)= -ve(:,1,2);
ne(:,1,3)= ve(:,2,3); ne(:,2,3)= -ve(:,1,3);
clear ve;
% for boundary edges e, neighbor(t,e) = t. So the difference is zero.     
edgeJump = dot((Du-Du(neighbor(:,1),:)),ne(:,:,1),2).^2 ...
         + dot((Du-Du(neighbor(:,2),:)),ne(:,:,2),2).^2 ...
         + dot((Du-Du(neighbor(:,3),:)),ne(:,:,3),2).^2;

%% Modification for Neumman boundary edges     
if (nargin == 5) && isfield(pde,'g_N')
    for k = 1:3
        idx = (bdFlag(:,k) == 2);	
        ne(idx,:,k) = ne(idx,:,k)./...
            repmat(sqrt(ne(idx,1,k).^2 + ne(idx,2,k).^2),1,2);  % normalize the normal vector ne 
        mid = (node(elem(idx,mod(k,3)+1),:) + node(elem(idx,mod(k+1,3)+1),:))/2;
        edgeJump(idx) = (pde.g_N(mid) - dot(Du(idx,:),ne(idx,:,k),2)).^2.*area(idx);        
    end
end         

%% Elementwise residual
elemResidual = zeros(size(elem,1),1);
if isfield(pde,'f') && isnumeric(pde.f) && (pde.f==0)
    pde.f = [];
end
if isfield(pde,'f') && ~isempty(pde.f)
    [lambda,weight] = quadpts(3);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        elemResidual = elemResidual + weight(p)*fp.^2;
    end
    elemResidual = elemResidual.*(area.^2);
end

%% Residual type error estimator
eta = (elemResidual + edgeJump).^(1/2);
% TODO an example with Neumann boundary condition
