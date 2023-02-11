function errK = getErrInf3D(uh, fun, fem, dind)

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

% Last Modified: 08/07/2020 by Xu Zhang 
%% 0. Initializaiton
if strcmp(fem.type,'P1')||strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR')
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem1.type,'P2')||strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end
%% 1. Error on noninterface elements
dof = fem.ldof; 
gx = fem.gx; gy = fem.gy; gz = fem.gz; 
uhK = uh(fem.t); 
tu = feval(fun, gx, gy, gz);
uK = 0; 
for i = 1:dof
    BAS = fem.bas(:,:,i);
    uK = uK + uhK(:,i).*feEvalBas(BAS, gx, gy, gz, dind);
end
errK = norm(tu-uK,'inf');
