function [err,errK] = getErrIFE3D(uh, pde, mesh, fem, femI, eNorm)

%% USAGE: calculate L^2 or semi-H^1 errors of FE/DG solution. 
%
% INPUTS:
% uh   --- a vector contains the FE solution or interpolation
% fun  --- exact function (u, Dx(u), Dy(u), or Dz(u))  
% mesh --- a struct data contains mesh information.
% fem  --- global degree of freedom for test function space
% type --- which norm 'L2','H1x','H1y','H1z'
%
% OUTPUTS:
% error --- error over the domain
%
% Last Modified: 09/18/2019 by Xu Zhang
% Last Modified: 07/02/2020 by Xu Zhang

%% 0. Initializaiton
if strcmp(fem.type,'P1') || strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR') 
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem.type,'P2') || strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end

if strcmp(eNorm,'L2')
    dind = [0,0,0];
    fun = pde.exactu;
    coefun = @(x,y,z) ones(size(x));
elseif strcmp(eNorm,'H1x') || strcmp(eNorm,'Ex')
    dind = [1,0,0];
    fun = pde.Dxu;
    coefun = @(x,y,z) ones(size(x));
elseif strcmp(eNorm,'H1y') || strcmp(eNorm,'Ey')
    dind = [0,1,0];
    fun = pde.Dyu;
    coefun = @(x,y,z) ones(size(x));
elseif strcmp(eNorm,'H1z') || strcmp(eNorm,'Ez')
    dind = [0,0,1];
    fun = pde.Dzu;
    coefun = @(x,y,z) ones(size(x));
end

if strcmp(eNorm,'Ex') || strcmp(eNorm,'Ey') || strcmp(eNorm,'Ez')
    coefun = pde.A;
end

%% 1. Error on non-interface elements
dof = fem.ldof;  
ntID = find(mesh.tLoc > 0); 
AN = fem.area(ntID); gw = fem.gw;
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); 

uhK = uh(fem.t); uhKn = uhK(ntID,:); 
tuN = feval(fun, gxN, gyN, gzN);
coefN = feval(coefun, gxN, gyN, gzN);
uKN = 0; errK = zeros(size(mesh.t,1),1);
for i = 1:dof
    BAS = fem.bas(ntID,:,i);
    uKN = uKN + uhKn(:,i).*feEvalBas(BAS, gxN, gyN, gzN, dind);
end
errKN = sqrt(AN.*sum(coefN.*(tuN-uKN).^2.*gw',2));
errK(ntID,1) = errKN;

%% 2. Error on interface elements 
A = femI.area; gw = femI.gw; gx = femI.gx; gy = femI.gy; gz = femI.gz; 
itID = find(mesh.tLoc < 0); ntI = length(itID);
uhK = uh(fem.t); uhKitmp = uhK(itID,:); uhKi = zeros(length(femI.t),4);
id = 0;
for i = 1:ntI
    tmp = uhKitmp(i,femI.locID(i,:));
    uhKi(id+1:id+femI.quadT(i),:) = repmat(tmp,femI.quadT(i),1);
    id = id + femI.quadT(i);
end
tuI = feval(fun, gx, gy, gz);
coefI = feval(coefun, gx, gy, gz);

uKI = 0; 
for i = 1:dof
    BAS = femI.bas(:,:,i);
    uKI = uKI + uhKi(:,i).*feEvalBas(BAS, gx, gy, gz, dind);
end
errKitmp = sqrt(A.*sum(coefI.*(tuI-uKI).^2.*gw',2));
errKI = zeros(ntI,1);

rID = 0;
for i = 1:ntI
    ntti = femI.quadT(i);
    errKI(i) = norm(errKitmp(rID+1:rID+ntti));
    rID = rID + ntti;
end
errK(itID) = errKI;
err = norm([errKN;errKI]);
