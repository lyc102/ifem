function err = getL2error3Q1(node,elem,exact_u,uh,option) 
%% GETL2ERROR3Q1 L2 norm of arroximation error on hex mesh.
%
%  err = getL2error(node,elem,@uexact,uh) computes the L2 norm of the error
%  between the exact solution uexact and a finite element approximation uh
%  on a hex mesh described by node and elem.
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Set up
if ~exist('option','var'), option = []; end
[N,dim] = size(node);
NT = size(elem,1);
if isfield('option','quadorder')
    quadOrder = option.quadorder; 
else
    quadOrder = 2;
end

%% Mesh size and Volume
hx = node(elem(:,2),1) - node(elem(:,1),1);
hy = node(elem(:,4),2) - node(elem(:,1),2);
hz = node(elem(:,5),3) - node(elem(:,1),3);
volume = abs(hx.*hy.*hz);

%%
[GaussPt, weight] = quadptscube(quadOrder);
lambdax = GaussPt(:,1);
lambday = GaussPt(:,2);
lambdaz = GaussPt(:,3);

nQuad = size(GaussPt,1);
err = zeros(NT,1);
for p = 1:nQuad
    if isfield(option,'isoparametric') && option.isoparametric
        [phi, tempvar, J] = hexbasis(node,elem,GaussPt(p,:));
    else
        phi2d = kron([lambday(p) 1-lambday(p)],[lambdax(p); 1-lambdax(p)]);
        phi = transpose(kron([lambdaz(p) 1-lambdaz(p)],phi2d([1 2 4 3])));
        J = volume;
    end    
%     [phi, tempvar, J] = hexbasis(node,elem, pts(p,:));
    pxyz = zeros(NT, dim);
    for i = 1:dim
        xi = node(:,i);
        pxyz(:,i) = xi(elem)*phi(:);
    end
    u = uh(elem)*phi(:);
    err = err + weight(p)*(u - exact_u(pxyz)).^2.*J;  
end
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));