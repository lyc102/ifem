function rhs = globRHSIFE3D(fun, mesh, fem, femI, dind)

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

%% 0. Initializaiton
if strcmp(fem.type,'P1') || strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR') 
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem.type,'P2') || strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end

%% 1. RHS on non-interface elements
dof = fem.ldof;  nloc = dof;
ntID = find(mesh.tLoc > 0); ntN = length(ntID);
AN = fem.area(ntID); gw = fem.gw;
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); 
X = zeros(nloc*ntN, 1);

fN = feval(fun,gxN,gyN,gzN);
ind = 0;
I = reshape(fem.t(ntID,:),nloc*ntN,1);
for i = 1:dof
    ibas = feEvalBas(fem.bas(ntID,:,i), gxN, gyN, gzN, dind);
    X(ind+1:ind+ntN) = AN.*sum((ibas.*fN).*gw',2);
    ind = ind + ntN;
end
rhsN = sparse(I,1,X,length(fem.p),1);

%% 2. RHS on interface elements
AI = femI.area; gw = femI.gw; gxI = femI.gx; gyI = femI.gy; gzI = femI.gz; 
ntI = size(femI.t,1); % not number of interface element, but quadrature element
X = zeros(nloc*ntI, 1);

fI = feval(fun,gxI,gyI,gzI);
ind = 0;
I = reshape(femI.t,nloc*ntI,1);
for i = 1:dof
    ibas = feEvalBas(femI.bas(:,:,i), gxI, gyI, gzI, dind);
    X(ind+1:ind+ntI) = AI.*sum((ibas.*fI).*gw',2);
    ind = ind + ntI;
end
rhsI = sparse(I,1,X,length(fem.p),1);
rhs = rhsN + rhsI;