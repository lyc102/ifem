function eta = estimateresidual3(node,elem,u,pde,bdFlag)
%% ESTIMATERESIDUAL3 residual type error estimator in 3-D.
%
% eta = estimateresidual3(node,elem,u,@f) computes the residual type a
% posterior error estimator for the Poisson equation with Dirichlet
% boundary condition or homogenous Neumann boundary condition g=0 in the
% form
%
% $$\eta _T^2 = volume(T)\sum _{F\in T}[\nabla u\cdot n_F]^2 + volume(T)^{5/3}\sum _i w_if_i^2,$$
%
% which is an approximation of 
%
% $$\eta _T^2 = \sum _{F\in T}\int _F h_F[\nabla u\cdot n_E]^2 dS + \int _T h^2f^2 dx.$$
% 
% Here [\nabla u\cdot n_F] denotes the jump of the flux accross the
% interior face F and n_F is the normal vector of F. In 3-D, Dlambda is an
% inward normal vector and n_F can be computed from normalization of
% Dlambda. In 3-D, h_F|area(F)| is approximated by volume(T).
%
% eta = estimateresidual3(node,elem,u,@f,pde) computes the estimator with
% more general setting. 
%
% On Neumann boundary faces, the jump as is modfied to g - du/dn and the
% contribution from the Neumann face is volume(T)*|g - du/dn|^2. 
%
% For general diffusion equation div(K grad u) = f, the flux is K grad u.
%
% PDE is a structure array that records the information of Neumann
% boundary condition and diffusion coefficent.
%
%  pde.g     : Neumann boundary condition data
%  pde.bdFace: boundary faces
%  pde.d     : diffusion coeffcients
%
% Example
%
% See also estimateresidual, estimaterecovery3
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.
 
%% Compute the flux from the finite element function u
[Du,volume,Dlambda] = gradu3(node,elem,u);
if exist('pde','var') && isfield(pde,'d') && ~isempty(pde.d)
    if isreal(pde.d)
        K = pde.d;                      % d is an array
    else                            % d is a function
        center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
                  node(elem(:,3),:) + node(elem(:,4),:))/4;
        K = pde.d(center);              
    end
    Du = [K.*Du(:,1),K.*Du(:,2),K.*Du(:,3)];
end
                             
%% Jump of nf flux
% data structure
T = auxstructure3(elem);
neighbor = T.neighbor;
clear T
% nf vector
nf(:,:,1) = -Dlambda(:,:,1)./repmat(sqrt(Dlambda(:,1,1).^2 + ...
                                 Dlambda(:,2,1).^2+Dlambda(:,3,1).^2),1,3);
nf(:,:,2) = -Dlambda(:,:,2)./repmat(sqrt(Dlambda(:,1,2).^2 + ...
                                 Dlambda(:,2,2).^2+Dlambda(:,3,2).^2),1,3);
nf(:,:,3) = -Dlambda(:,:,3)./repmat(sqrt(Dlambda(:,1,3).^2 + ...
                                 Dlambda(:,2,3).^2+Dlambda(:,3,3).^2),1,3);
nf(:,:,4) = -Dlambda(:,:,4)./repmat(sqrt(Dlambda(:,1,4).^2 + ...
                                 Dlambda(:,2,4).^2+Dlambda(:,3,4).^2),1,3);
% for boundary edges e, neighbor(t,e) = t. So the difference is zero.     
faceJump = dot((Du-Du(neighbor(:,1),:)),nf(:,:,1),2).^2 ...
         + dot((Du-Du(neighbor(:,2),:)),nf(:,:,2),2).^2 ...
         + dot((Du-Du(neighbor(:,3),:)),nf(:,:,3),2).^2 ...
         + dot((Du-Du(neighbor(:,4),:)),nf(:,:,4),2).^2;
faceJump = volume.*faceJump;

%% Modification for Neumman boundary edges     
if exist('pde','var') && isfield(pde,'gN') && ~isempty(pde.gN)
    for k = 1:4
        idx = (bdFlag(:,k) == 2);	        
        mid = (node(elem(idx,1),:) + node(elem(idx,2),:) + node(elem(idx,3),:) ...
             + node(elem(idx,4),:) - node(elem(idx,k),:))/3;
        faceJump(idx) = (pde.gN(mid) - dot(Du(idx),nf(idx,:,k),2)).^2.*volume(idx);        
    end
end         

%% Compute element residual using 4-points rule
elemResidual = zeros(size(elem,1),1);
if isfield(pde,'f') && isnumeric(pde.f) && (pde.f==0)
    pde.f = [];
end
if isfield(pde,'f') && ~isempty(pde.f)
    [lambda,weight] = quadpts3(3);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxyz = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:) ...
            + lambda(p,4)*node(elem(:,4),:);
        fp = pde.f(pxyz);
        elemResidual = elemResidual + weight(p)*fp.^2;
    end
    elemResidual = elemResidual.*volume.^(5/3);
end

%% Compute residual type error estimator
eta = (elemResidual + faceJump).^(1/2);