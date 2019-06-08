function err = getL2errorQ1(node,elem, exact_u, uh) 
%% GETL2ERRORQ1 L2 norm of arroximation error on quad mesh.
%
%  err = getL2error(node,elem,@uexact,uh) computes the L2 norm of the error
%  between the exact solution uexact and a finite element approximation uh
%  on a quad mesh described by node and elem.
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[N,  n] = size(node);
NT = size(elem,1);

[pts, weight] = quadquadpts(3);
nQuad = size(pts,1);
err = zeros(NT,1);
for p = 1:nQuad
    [phi, tempvar, J] = quadbasis(node,elem, pts(p,:));
    pxy = zeros(NT, n);
    for i = 1:n
        xi = node(:,i);
        pxy(:,i) = xi(elem)*phi;
    end
    u = uh(elem)*phi;
    err = err + weight(p)*(u - exact_u(pxy)).^2.*J;  
end

err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));