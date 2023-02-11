function S = globMatrixNedVIFE3D(fun,MatInd,mesh,femI,fem1,fem2)

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
feEvalBas1 = @EvalNed1Bas3D;
feEvalBas2 = @EvalNed1Bas3D;

IfeEvalBas1 = @evalNed1IFEBas3D;
IfeEvalBas2 = @evalNed1IFEBas3D;

dof1 = fem1.ldof; dof2 = fem2.ldof; nloc = dof1*dof2; 
ntID = find(mesh.tLoc > 0); ntN = length(ntID);
t_e_orit = fem1.t_e_orit;

%% 1. Matrix on noninterface elements
AN = fem1.area(ntID); 
gxN = fem1.gx(ntID,:); gyN = fem1.gy(ntID,:); gzN = fem1.gz(ntID,:); gw = fem1.gw;
XN = zeros(nloc*ntN, 1);

coefN = feval(fun,gxN,gyN,gzN);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);

% nt = length(mesh.t);
% A = fem1.area; gx = fem1.gx; gy = fem1.gy; gz = fem1.gz; gw = fem1.gw;
% X = zeros(nloc*nt, 1);
% coef = feval(fun,gx,gy,gz);
% Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
% Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);
for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem1.bas, ntID, gxN, gyN, gzN, i, MatInd, 1).*t_e_orit(ntID,i);
    Ibasy{i} = feEvalBas1(fem1.bas, ntID, gxN, gyN, gzN, i, MatInd, 2).*t_e_orit(ntID,i);
    Ibasz{i} = feEvalBas1(fem1.bas, ntID, gxN, gyN, gzN, i, MatInd, 3).*t_e_orit(ntID,i);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem2.bas, ntID, gxN, gyN, gzN, j, MatInd, 1).*t_e_orit(ntID,j);
    Jbasy{j} = feEvalBas2(fem2.bas, ntID, gxN, gyN, gzN, j, MatInd, 2).*t_e_orit(ntID,j);
    Jbasz{j} = feEvalBas2(fem2.bas, ntID, gxN, gyN, gzN, j, MatInd, 3).*t_e_orit(ntID,j);
end

IN = reshape(repmat(femI.g2ldof(ntID,1:6),6,1),nloc*ntN,1);
JN = repmat(reshape(femI.g2ldof(ntID,1:6),dof2*ntN,1),6,1);
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
SN = sparse(IN(ID),JN(ID),XN(ID),size(femI.gdof,1),size(femI.gdof,1));

%% 2. Matrix on interface elements

% II = femI.II; JI = femI.JI; XI = femI.XI;
% SI = sparse(II,JI,XI,size(femI.g2ldof,1),size(femI.g2ldof,1));

S = SN;% + SI;