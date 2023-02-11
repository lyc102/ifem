function rhs = globNedFitRHS3D(fm,fp, mesh, fem, dind, vind)

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

feEvalBas = @EvalNed1Bas3D;

f = feval(fm,gx,gy,gz);
tId2 = (mesh.tLoc == 2);
f(tId2,:) = feval(fp,gx(tId2,:),gy(tId2,:),gz(tId2,:));
ind = 0;
I = reshape(fem.g2ldof,nloc*nt,1);
for i = 1:dof
    ibas = feEvalBas(fem.bas, ':', gx, gy, gz, i, dind, vind);
    X(ind+1:ind+nt) = A.*sum((ibas.*f).*gw',2).*fem.t_e_orit(:,i);
    ind = ind + nt;
end
rhs = sparse(I,1,X,length(fem.gdof),1);