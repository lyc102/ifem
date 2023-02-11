function S = globMatrixVIFE3DMass(fun,mesh,fem1,fem2)

%% USAGE: generate stiffness global matrix on a 3D mesh 
%
% INPUTS:
% fun --- coefficient function
% mesh --- a struct data contains very rich mesh information.
% fem1 --- global DoF for test function space
% fem2 --- global DoF for trial function space
%
% OUTPUTS:
% [IN JN XN] --- triplets of the sparse matrix from regular elements. 
% [II JI XI] --- triplets of the sparse matrix from interface elements. 

% Last Modified: 08/07/2020 by Xu Zhang 

%% 0. Initializaiton
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
Ibas = cell(dof1,1); 
Jbas = cell(dof2,1);
for i = 1:dof1
    Ibas{i} = feEvalBas1(fem1.bas(ntID,:,i), gxN, gyN, gzN, [0,0,0]);
end
for j = 1:dof2
    Jbas{j} = feEvalBas2(fem2.bas(ntID,:,j), gxN, gyN, gzN, [0,0,0]); 
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
ID = find(XN~=0); 
Ndof = size([mesh.p;mesh.eIntP],1);
S = sparse(IN(ID),JN(ID),XN(ID),Ndof,Ndof);