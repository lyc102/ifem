function [IN,JN,XN,II,JI,XI] = globMatrixIFE3D(fun,mesh,femI,fem1,d1,fem2,d2)

%% USAGE: generate global matrix on a 3D mesh
%
% INPUTS:
% fun --- coefficient function
% mesh --- a struct data contains very rich mesh information.
% fem1 --- global DoF for test function space
% fem2 --- global DoF for trial function space
% d1 --- derivative info for test function
% d2 --- derivative info for trial function
%            d = [0,0,0]: function value
%            d = [1,0,0]: Dx value
%            d = [0,1,0]: Dy value
%            d = [0,0,1]: Dz value
%
% OUTPUTS:
% [IN JN XN] --- triplets of the sparse matrix from regular elements. 
% [II JI XI] --- triplets of the sparse matrix from interface elements. 

% Last Modified: 08/07/2020 by Xu Zhang 

%% 0. Initializaiton
if nargin == 5
    fem2 = fem1; d2 = d1;
end

if strcmp(fem1.type,'P1')||strcmp(fem1.type,'DGP1')||strcmp(fem1.type,'CR')
    feEvalBas1 = @evalP1Bas3D;
elseif strcmp(fem1.type,'P2')||strcmp(fem1.type,'DGP2')
    feEvalBas1 = @evalP2Bas3D;
end

if strcmp(fem2.type,'P1')||strcmp(fem2.type,'DGP1')||strcmp(fem2.type,'CR')
    feEvalBas2 = @evalP1Bas3D;
elseif strcmp(fem2.type,'P2')||strcmp(fem2.type,'DGP2')
    feEvalBas2 = @evalP2Bas3D;
end

%% 1. Matrix on noninterface elements
dof1 = fem1.ldof; dof2 = fem2.ldof; nloc = dof1*dof2; 
ntID = find(mesh.tLoc > 0); ntN = length(ntID);
AN = fem1.area(ntID); 
gxN = fem1.gx(ntID,:); gyN = fem1.gy(ntID,:); gzN = fem1.gz(ntID,:); gw = fem1.gw;
XN = zeros(nloc*ntN, 1);

coefN = feval(fun,gxN,gyN,gzN);
Ibas = cell(dof1,1); Jbas = cell(dof2,1);
for i = 1:dof1
    Ibas{i} = feEvalBas1(fem1.bas(ntID,:,i), gxN, gyN, gzN, d1);
end
for j = 1:dof2
    Jbas{j} = feEvalBas2(fem2.bas(ntID,:,j), gxN, gyN, gzN, d2);
end

IN = reshape(repmat(fem1.t(ntID,:),4,1),nloc*ntN,1);
JN = repmat(reshape(fem2.t(ntID,:),dof2*ntN,1),4,1);
ind = 0;
for i = 1:dof1
    for j = 1:dof2
        XN(ind+1:ind+ntN) = AN.*sum(((Ibas{i}.*(coefN.*Jbas{j})).*gw'),2);
        ind = ind + ntN;
    end
end

%% 2. Matrix on interface elements
AI = femI.area; gw = femI.gw; gxI = femI.gx; gyI = femI.gy; gzI = femI.gz; 
ntI = size(femI.t,1); % not number of interface element, but quadrature element
XI = zeros(nloc*ntI, 1);

coefI = feval(fun,gxI,gyI,gzI);
Ibas = cell(dof1,1); Jbas = cell(dof2,1);
for i = 1:dof1
    Ibas{i} = feEvalBas1(femI.bas(:,:,i), gxI, gyI, gzI, d1);
end
for j = 1:dof2
    Jbas{j} = feEvalBas2(femI.bas(:,:,j), gxI, gyI, gzI, d2);
end

II = reshape(repmat(femI.t,4,1),nloc*ntI,1);
JI = repmat(reshape(femI.t,dof2*ntI,1),4,1);
ind = 0; 
for i = 1:dof1
    for j = 1:dof2
        XI(ind+1:ind+ntI) = AI.*sum(((Ibas{i}.*(coefI.*Jbas{j})).*gw'),2);
        ind = ind + ntI;
    end
end