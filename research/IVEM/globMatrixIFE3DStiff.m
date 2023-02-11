function S = globMatrixIFE3DStiff(fun,mesh,femI,fem1,fem2)

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
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);
for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem1.bas(ntID,:,i), gxN, gyN, gzN, [1,0,0]);
    Ibasy{i} = feEvalBas1(fem1.bas(ntID,:,i), gxN, gyN, gzN, [0,1,0]);
    Ibasz{i} = feEvalBas1(fem1.bas(ntID,:,i), gxN, gyN, gzN, [0,0,1]);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem2.bas(ntID,:,j), gxN, gyN, gzN, [1,0,0]); 
    Jbasy{j} = feEvalBas2(fem2.bas(ntID,:,j), gxN, gyN, gzN, [0,1,0]);
    Jbasz{j} = feEvalBas2(fem2.bas(ntID,:,j), gxN, gyN, gzN, [0,0,1]);
end

IN = reshape(repmat(fem1.t(ntID,:),4,1),nloc*ntN,1);
JN = repmat(reshape(fem2.t(ntID,:),dof2*ntN,1),4,1);
ind = 0;
for i = 1:dof1
    for j = 1:dof2
        XN(ind+1:ind+ntN) = AN.*(sum(((Ibasx{i}.*(coefN.*Jbasx{j})).*gw'),2) + ...
            sum(((Ibasy{i}.*(coefN.*Jbasy{j})).*gw'),2) + ...
            sum(((Ibasz{i}.*(coefN.*Jbasz{j})).*gw'),2));
        ind = ind + ntN;
    end
end
ID = find(XN~=0); 
SN = sparse(IN(ID),JN(ID),XN(ID),size(fem1.p,1),size(fem2.p,1));

%% 2. Matrix on interface elements
AI = femI.area; gw = femI.gw; gxI = femI.gx; gyI = femI.gy; gzI = femI.gz; 
ntI = size(femI.t,1); % not number of interface element, but quadrature element
XI = zeros(nloc*ntI, 1);

coefI = feval(fun,gxI,gyI,gzI);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);
for i = 1:dof1
    Ibasx{i} = feEvalBas1(femI.bas(:,:,i), gxI, gyI, gzI, [1,0,0]);
    Ibasy{i} = feEvalBas1(femI.bas(:,:,i), gxI, gyI, gzI, [0,1,0]);
    Ibasz{i} = feEvalBas1(femI.bas(:,:,i), gxI, gyI, gzI, [0,0,1]);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(femI.bas(:,:,j), gxI, gyI, gzI, [1,0,0]);
    Jbasy{j} = feEvalBas2(femI.bas(:,:,j), gxI, gyI, gzI, [0,1,0]);
    Jbasz{j} = feEvalBas2(femI.bas(:,:,j), gxI, gyI, gzI, [0,0,1]);
end

II = reshape(repmat(femI.t,4,1),nloc*ntI,1);
JI = repmat(reshape(femI.t,dof2*ntI,1),4,1);
ind = 0; 
for i = 1:dof1
    for j = 1:dof2
        XI(ind+1:ind+ntI) = AI.*(sum(((Ibasx{i}.*(coefI.*Jbasx{j})).*gw'),2) + ...
            sum(((Ibasy{i}.*(coefI.*Jbasy{j})).*gw'),2) + ...
            sum(((Ibasz{i}.*(coefI.*Jbasz{j})).*gw'),2));
        ind = ind + ntI;
    end
end
SI = sparse(II,JI,XI,size(fem1.p,1),size(fem2.p,1));
S = SN + SI;