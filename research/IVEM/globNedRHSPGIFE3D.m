function rhs = globNedRHSPGIFE3D(fun, mesh, fem, femI, dind, vind)

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
dof = size(fem.g2ldof,2);  nloc = dof;
ntID = find(mesh.tLoc > 0); ntN = length(ntID);
AN = fem.area(ntID); gw = fem.gw;
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); 
X = zeros(nloc*ntN, 1);

feEvalBas = @EvalNed1Bas3D;

fN = feval(fun,gxN,gyN,gzN);
ind = 0;
I = reshape(fem.g2ldof(ntID,:),nloc*ntN,1);
for i = 1:dof
    ibas = feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, vind);
    X(ind+1:ind+ntN) = AN.*sum((ibas.*fN).*gw',2).*fem.t_e_orit(ntID,i);
    ind = ind + ntN;
end
rhsN = sparse(I,1,X,length(fem.gdof),1);

%% 2. RHS on interface elements
AI = femI.area; gw = femI.gw; gxI = femI.gx; gyI = femI.gy; gzI = femI.gz; 
ntI = size(femI.g2ldof,1); % not number of interface element, but quadrature element
X = zeros(nloc*ntI, 1);
tIntfID = femI.tIntfID;

fI = feval(fun,gxI,gyI,gzI);
ind = 0;
I = reshape(femI.g2ldof,nloc*ntI,1);
for i = 1:dof
    ibas = feEvalBas(fem.bas,tIntfID, gxI, gyI, gzI,i, dind, vind).*fem.t_e_orit(tIntfID,i);
    X(ind+1:ind+ntI) = AI.*sum((ibas.*fI).*gw',2);
    ind = ind + ntI;
end
rhsI = sparse(I,1,X,size(fem.gdof,1),1);
rhs = rhsN + rhsI;
