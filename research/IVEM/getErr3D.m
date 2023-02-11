function [err,errK] = getErr3D(uh, fun, fem, eNorm)

%% USAGE: calculate L^2 or semi-H^1 errors of FE/DG solution. 
%
% INPUTS:
% uh   --- a vector contains the FE solution or interpolation
% fun  --- exact function (u, Dx(u), Dy(u), or Dz(u))  
% fem  --- global degree of freedom for test function space
% type --- which norm 'L2','H1x','H1y','H1z'
%
% OUTPUTS:
% error --- error over the domain
%
% Last Modified: 08/07/2020 by Xu Zhang

%% 1. Error on non-interface elements 
dof = fem.ldof; 
A = fem.area; gw = fem.gw; gx = fem.gx; gy = fem.gy; gz = fem.gz; 
tu = feval(fun, gx, gy, gz);
uhK = uh(fem.t); 
if strcmp(fem.type,'P1') || strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR') 
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem.type,'P2') || strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end

if strcmp(eNorm,'L2')
    dind = [0,0,0];
elseif strcmp(eNorm,'H1x')
    dind = [1,0,0];
elseif strcmp(eNorm,'H1y')
    dind = [0,1,0];
elseif strcmp(eNorm,'H1z')
    dind = [0,0,1];
end
uK = 0;
for i = 1:dof
    BAS = fem.bas(:,:,i);
    uK = uK + uhK(:,i).*feEvalBas(BAS, gx, gy, gz, dind);
end
errK = sqrt(A.*sum((tu-uK).^2.*gw',2));
err = norm(errK);
