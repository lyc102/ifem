function [err,errK1,errK2,errK3] = getCurlErr3D(uh, pde, fem, eNorm)

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

if strcmp(eNorm,'L2')
    dind = 0;
    fun1 = pde.exactu1; fun2 = pde.exactu2; fun3 = pde.exactu3;
elseif strcmp(eNorm,'Curl')
    dind = 1;
    fun1 = pde.Dxu; fun2 = pde.Dyu; fun3 = pde.Dzu;
end
tu1 = feval(fun1, gx, gy, gz);
tu2 = feval(fun2, gx, gy, gz);
tu3 = feval(fun3, gx, gy, gz);
uhK = uh(fem.g2ldof); 
feEvalBas = @EvalNed1Bas3D;


uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:dof
    uK1 = uK1 + uhK(:,i).*feEvalBas(fem.bas, ':', gx, gy, gz, i, dind, 1).*fem.t_e_orit(:,i);
    uK2 = uK2 + uhK(:,i).*feEvalBas(fem.bas, ':', gx, gy, gz, i, dind, 2).*fem.t_e_orit(:,i);
    uK3 = uK3 + uhK(:,i).*feEvalBas(fem.bas, ':', gx, gy, gz, i, dind, 3).*fem.t_e_orit(:,i);
end
errK1 = sqrt(A.*sum((tu1-uK1).^2.*gw',2));
errK2 = sqrt(A.*sum((tu2-uK2).^2.*gw',2));
errK3 = sqrt(A.*sum((tu3-uK3).^2.*gw',2));


err = sum((errK1.^2+errK2.^2+errK3.^2)).^(1/2);
