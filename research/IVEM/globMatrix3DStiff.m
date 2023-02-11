function S = globMatrix3DStiff(fun,mesh,fem1,fem2)

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

dof1 = fem1.ldof; dof2 = fem2.ldof; nloc = dof1*dof2; 
nt = length(mesh.t);
A = fem1.area; gx = fem1.gx; gy = fem1.gy; gz = fem1.gz; gw = fem1.gw;
X = zeros(nloc*nt, 1);

coef = feval(fun,gx,gy,gz);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);
for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem1.bas(:,:,i), gx, gy, gz, [1,0,0]);
    Ibasy{i} = feEvalBas1(fem1.bas(:,:,i), gx, gy, gz, [0,1,0]);
    Ibasz{i} = feEvalBas1(fem1.bas(:,:,i), gx, gy, gz, [0,0,1]);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem2.bas(:,:,j), gx, gy, gz, [1,0,0]); 
    Jbasy{j} = feEvalBas2(fem2.bas(:,:,j), gx, gy, gz, [0,1,0]);
    Jbasz{j} = feEvalBas2(fem2.bas(:,:,j), gx, gy, gz, [0,0,1]);
end

I = reshape(repmat(fem1.t,4,1),nloc*nt,1);
J = repmat(reshape(fem2.t,dof2*nt,1),4,1);
ind = 0;
for i = 1:dof1
    for j = 1:dof2
        X(ind+1:ind+nt) = A.*(sum(((Ibasx{i}.*(coef.*Jbasx{j})).*gw'),2) + ...
            sum(((Ibasy{i}.*(coef.*Jbasy{j})).*gw'),2) + ...
            sum(((Ibasz{i}.*(coef.*Jbasz{j})).*gw'),2));
        ind = ind + nt;
    end
end
ID = find(X~=0); 
S = sparse(I(ID),J(ID),X(ID),size(fem1.p,1),size(fem2.p,1));
