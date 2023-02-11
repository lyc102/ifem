function [err,errK1,errK2,errK3] = getCurlErrIFE3D(uh, pde, mesh, fem, femI, eNorm)

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

feEvalBas = @EvalNed1Bas3D;
IfeEvalBas = @evalNed1IFEBas3D;

%% 1. Error on non-interface elements 
dof = fem.ldof;  
ntID = find(mesh.tLoc>0); 
AN = fem.area(ntID); gw = fem.gw;
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:);
errK1 = zeros(size(mesh.t,1),1); errK2 = errK1; errK3 = errK1;

if strcmp(eNorm,'L2')
    dind = 0;
    fun1 = pde.exactu1; fun2 = pde.exactu2; fun3 = pde.exactu3;
elseif strcmp(eNorm,'Curl')
    dind = 1;
    fun1 = pde.Dxu; fun2 = pde.Dyu; fun3 = pde.Dzu;
end
tu1 = feval(fun1, gxN, gyN, gzN);
tu2 = feval(fun2, gxN, gyN, gzN);
tu3 = feval(fun3, gxN, gyN, gzN);
uhK = uh(fem.g2ldof); 


uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:dof
    uK1 = uK1 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, 1).*fem.t_e_orit(ntID,i);
    uK2 = uK2 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, 2).*fem.t_e_orit(ntID,i);
    uK3 = uK3 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, 3).*fem.t_e_orit(ntID,i);
end
errK1(ntID) = sqrt(AN.*sum((tu1-uK1).^2.*gw',2));
errK2(ntID) = sqrt(AN.*sum((tu2-uK2).^2.*gw',2));
errK3(ntID) = sqrt(AN.*sum((tu3-uK3).^2.*gw',2));

%err = sum((errK1.^2+errK2.^2+errK3.^2)).^(1/2);


%% 1. Error on interface elements 
A = femI.area; gw = femI.gw; gx = femI.gx; gy = femI.gy; gz = femI.gz; 
tIntfID = femI.tIntfID;
itID = find(mesh.tLoc < 0); ntI = length(itID);
uhK = uh(fem.g2ldof); uhKitmp = uhK(itID,:); uhKi = zeros(length(femI.g2ldof),6);
id = 0;
for i = 1:ntI
    tmp = uhKitmp(i,:);
    uhKi(id+1:id+femI.quadT(i),:) = repmat(tmp,femI.quadT(i),1);
    id = id + femI.quadT(i);
end
tu1 = feval(fun1, gx, gy, gz);
tu2 = feval(fun2, gx, gy, gz);
tu3 = feval(fun3, gx, gy, gz);

uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:dof
    BAS = femI.bas(:,:,i);
    uK1 = uK1 + uhKi(:,i).*IfeEvalBas(BAS, gx, gy, gz, dind, 1).*fem.t_e_orit(tIntfID,i);
    uK2 = uK2 + uhKi(:,i).*IfeEvalBas(BAS, gx, gy, gz, dind, 2).*fem.t_e_orit(tIntfID,i);
    uK3 = uK3 + uhKi(:,i).*IfeEvalBas(BAS, gx, gy, gz, dind, 3).*fem.t_e_orit(tIntfID,i);
end
errKitmp1 = sqrt(A.*sum((tu1-uK1).^2.*gw',2));
errKitmp2 = sqrt(A.*sum((tu2-uK2).^2.*gw',2));
errKitmp3 = sqrt(A.*sum((tu3-uK3).^2.*gw',2));

errKI1 = zeros(ntI,1); errKI2 = errKI1; errKI3 = errKI1;

rID = 0;
for i = 1:ntI
    ntti = femI.quadT(i);
    errKI1(i) = norm(errKitmp1(rID+1:rID+ntti));
    errKI2(i) = norm(errKitmp2(rID+1:rID+ntti));
    errKI3(i) = norm(errKitmp3(rID+1:rID+ntti));
    rID = rID + ntti;
end
errK1(itID) = errKI1;
errK2(itID) = errKI2;
errK3(itID) = errKI3;

err = sum((errK1.^2+errK2.^2+errK3.^2)).^(1/2);%*ntI^(1/4);




