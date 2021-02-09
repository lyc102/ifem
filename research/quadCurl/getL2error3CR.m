function [err, errK] = getL2error3CR(node,elem,uvecexact,uh,quadOrder)
%% GETL2ERROR3CR L2 norm of the approximation error for vectorial CR element
%    need Matlab ver > 2016b for the vectorization to work
%    
%    Reference: Error analysis of a decoupled finite element method for quad-curl problems
%               https://arxiv.org/abs/2102.03396
%
%    see alos: getL2error3
%  
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('quadOrder','var'); quadOrder = 3; end
elem2face = dof3face(elem);

%% 
NT = size(elem,1);
err = zeros(NT,1);
[lambda,weight] = quadpts3(quadOrder);
phi = 1-3*lambda;

nQuad = size(lambda,1);
for p = 1:nQuad
    uhp = uh(elem2face(:,1),:)*phi(p,1) + ...
        uh(elem2face(:,2),:)*phi(p,2) + ...
        uh(elem2face(:,3),:)*phi(p,3) + ...
        uh(elem2face(:,4),:)*phi(p,4);
    % quadrature points in the x-y coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
         + lambda(p,2)*node(elem(:,2),:) ...
         + lambda(p,3)*node(elem(:,3),:) ...
         + lambda(p,4)*node(elem(:,4),:);
    err = err + weight(p)*sum((uvecexact(pxyz) - uhp).^2,2);
end

%% 
d12 = node(elem(:,2),:)-node(elem(:,1),:);
d13 = node(elem(:,3),:)-node(elem(:,1),:);
d14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = abs(dot(mycross(d12,d13,2),d14,2)/6);
err(isnan(err)) = 0; % singular point is excluded
err = volume.*err;
if nargout > 1; errK = sqrt(err); end
err = sqrt(sum(err));

end