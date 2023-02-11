function [errKN,errKI] = getErrInfIFE3D(uh, fun, mesh, femI, fem, dind)

%% USAGE: generate u-uh on all Gaussian points 
%
% INPUTS:
% fun --- coefficient function
% mesh --- a struct data contains very rich mesh information.
% fem --- global DoF for test function space
% femI --- global DoF for test function space
% dind --- derivative info for test function
%            dind = [0,0,0]: function value
%            dind = [1,0,0]: Dx value
%            dind = [0,1,0]: Dy value
%            dind = [0,0,1]: Dz value
%
% OUTPUTS:
% matrix --- global (mass, stiffness, ...) matrix.

% Last Modified: 07/11/2020 by Xu Zhang 
%% 0. Initializaiton
if strcmp(fem.type,'P1')||strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR')
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem1.type,'P2')||strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end
%% 1. Error on noninterface elements
dof = fem.ldof; 
gx = fem.gx; gy = fem.gy; gz = fem.gz; 
ntID = find(mesh.tLoc > 0); 
gxN = gx(ntID,:); gyN = gy(ntID,:); gzN = gz(ntID,:); 
uhK = uh(fem.t); uhKN = uhK(ntID,:); 
tuN = feval(fun, gxN, gyN, gzN);
uKN = 0; 
for i = 1:dof
    BAS = fem.bas(ntID,:,i);
    uKN = uKN + uhKN(:,i).*feEvalBas(BAS, gxN, gyN, gzN, dind);
end
errKN = norm(tuN-uKN,'inf');

%% 2. Error on interface elements 
gxI = femI.gx; gyI = femI.gy; gzI = femI.gz; 
itID = find(mesh.tLoc < 0); ntI = length(itID);
uhK = uh(fem.t); uhKitmp = uhK(itID,:); uhKI = zeros(length(femI.t),4);
id = 0;
for i = 1:ntI
    tmp = uhKitmp(i,femI.locID(i,:));
    uhKI(id+1:id+femI.quadT(i),:) = repmat(tmp,femI.quadT(i),1);
    id = id + femI.quadT(i);
end
tuI = feval(fun,gxI,gyI,gzI);
uKI = 0; 
for i = 1:dof
    BAS = femI.bas(:,:,i);
    uKI = uKI + uhKI(:,i).*feEvalBas(BAS, gxI, gyI, gzI, dind);
end
errKI = norm(tuI-uKI,'inf');