function err = getH1errorQ1(node,elem,Du,uh,K,quadOrder)
%% GETH1ERRORQ1 H1 norm of the approximation error on quad mesh.
%
%  err = getH1errorQ1(node,elem,@Du,uh) computes the H1 norm of the
%  error between the exact solution Du and finite element approximation
%  uh on a quad mesh described by node and elem. 
%
%  The input parameter Du is a function handle uh is a column array with
%  length N representing a piecewise linear function on the mesh.
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
% Modified by Long Chen. The original computation is not quite right. The
% diffusion coefficient cannot be computed separately.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[N,  dim] = size(node);
[NT, nv] = size(elem);
if ~exist('quadOrder','var')
    quadOrder = 3; 
end

%% compute gradient of basis.
[pts, weight] = quadquadpts(quadOrder);
nQuad = size(pts,1);
err = zeros(size(elem,1),1);
for p = 1:nQuad
    [phi, Dphi, J] = quadbasis(node,elem, pts(p,:));
    pxy = zeros(NT, dim);
    for i = 1:dim
        xi = node(:,i);
        pxy(:,i) = xi(elem)*phi;
    end
    Duh = zeros(NT,dim);
    for i = 1:dim
        Dphi_i(:,1:nv) = Dphi(:,i,:);
        Duh(:,i) = sum(uh(elem).*Dphi_i,2);
    end
    if exist('K','var') && ~isempty(K) && ~isnumeric(K) % K is a function
        err = err + weight(p)*K(pxy).*sum((Du(pxy)-Duh).^2,2).*J;
    else
        err = err + weight(p)*sum((Du(pxy)-Duh).^2,2).*J;        
    end
end
if exist('K','var') && ~isempty(K) && isnumeric(K) && size(K,1) == NT
    err = K.*err;    % K is piecewise constant
end
%err = err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));
end
%% TODO: add more help information. 