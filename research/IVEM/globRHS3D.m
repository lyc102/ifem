function rhs = globRHS3D(fun, mesh, fem, dind)

%% USAGE: generate global load vector on a tetrahedral mesh
%
% INPUTS:
% fun --- the load function from PDE (pde.f)
% mesh --- a struct data contains mesh information.
% fem --- global DoF for test function space
% dind --- derivative info for test function
%            d = [0,0,0]: function value
%            d = [1,0,0]: Dx value
%            d = [0,1,0]: Dy value
%            d = [0,0,1]: Dz value
%
% OUTPUTS:
% rhs --- global rhs vector

% Last Modified: 08/07/2020 by Xu Zhang

%% 1. RHS on non-interface elements
dof = fem.ldof;  nloc = dof; nt = length(mesh.t);
A = fem.area; gw = fem.gw; gx = fem.gx; gy = fem.gy; gz = fem.gz;
X = zeros(nloc*nt, 1);

if strcmp(fem.type,'P1') || strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR')
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem.type,'P2') || strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end

f = feval(fun,gx,gy,gz);
ind = 0;
I = reshape(fem.t,nloc*nt,1);
for i = 1:dof
    ibas = feEvalBas(fem.bas(:,:,i), gx, gy, gz, dind);
    X(ind+1:ind+nt) = A.*sum((ibas.*f).*gw',2);
    ind = ind + nt;
end
rhs = sparse(I,1,X,length(fem.p),1);