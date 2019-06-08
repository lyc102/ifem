function err = getH1errorbd(node,elem,pde,uh,quadOrder)
%% GETH1ERRORBD H1 norm of the approximation error through boundary integral.
%
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('quadOrder','var')
    quadOrder = 3; 
end
Nu = size(uh,1);    N = size(node,1);   NT = size(elem,1);

[bdNode,bdEdge] = findboundary(elem);
err = zeros(NT,1);
[lambda,weight] = quadpts(quadOrder);

err = zeros(NT,1);

%% Element-wise integral
if isfield(pde,'f') && isnumeric(pde.f) && (pde.f==0)
    pde.f = [];
end
if isfield(pde,'f') && ~isempty(pde.f)
   [lambda,weight] = quadpts(quadOrder);
   phi = lambda; % linear bases 
   nQuad = size(lambda,1);
   for p = 1:nQuad
      uhp = uh(elem(:,1))*phi(p,1) + ...
            uh(elem(:,2))*phi(p,2) + ...
            uh(elem(:,3))*phi(p,3);
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
        
      fp = pde.f(pxy);
   end

end
err = area.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));