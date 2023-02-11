function rhs = globNedRHSVIFE3D(fun, mesh, fem, femI, dind, vind)

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
IfeEvalBas = @evalNed1IFEBas3D;

fN = feval(fun,gxN,gyN,gzN);
ind = 0;
I = reshape(femI.g2ldof(ntID,1:6),nloc*ntN,1);
for i = 1:dof
    ibas = feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, vind);
    X(ind+1:ind+ntN) = AN.*sum((ibas.*fN).*gw',2).*fem.t_e_orit(ntID,i);
    ind = ind + ntN;
end
rhsN = sparse(I,1,X,length(femI.gdof),1);

%% 2. RHS on interface elements
AI1 = femI.area1; gw = femI.gw; gxI1 = femI.gx1; gyI1 = femI.gy1; gzI1 = femI.gz1;
ntI1 = size(femI.tType1,1); % not number of interface element, but quadrature element
nloc1 = 12; X1 = zeros(nloc1*ntI1, 1);
fI1 = feval(fun,gxI1,gyI1,gzI1);
ind = 0;
I1 = reshape(femI.tType1,nloc1*ntI1,1);
for i = 1:nloc1
    ibas = femI.basUType1(:,vind,i);
    X1(ind+1:ind+ntI1) = AI1.*sum((ibas.*fI1).*gw',2);
    ind = ind + ntI1;
end
rhsI1 = sparse(I1,1,X1,size(femI.gdof,1),1);



AI2 = femI.area2; gw = femI.gw; gxI2 = femI.gx2; gyI2 = femI.gy2; gzI2 = femI.gz2;
ntI2 = size(femI.tType2,1); % not number of interface element, but quadrature element
nloc2 = 14; X2 = zeros(nloc2*ntI2, 1);
fI2 = feval(fun,gxI2,gyI2,gzI2);
ind = 0;
I2 = reshape(femI.tType2,nloc2*ntI2,1);
for i = 1:nloc2
    ibas = femI.basUType2(:,vind,i);
    X2(ind+1:ind+ntI2) = AI2.*sum((ibas.*fI2).*gw',2);
    ind = ind + ntI2;
end
rhsI2 = sparse(I2,1,X2,size(femI.gdof,1),1);


rhs = rhsN + rhsI1 + rhsI2;
