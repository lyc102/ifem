function S = globMatrixNed3D(fun,MatInd,mesh,fem1,fem2)

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

dof1 = fem1.ldof; dof2 = fem2.ldof; nloc = dof1*dof2; 
nt = length(mesh.t);
A = fem1.area; gx = fem1.gx; gy = fem1.gy; gz = fem1.gz; gw = fem1.gw;
X = zeros(nloc*nt, 1);

coef = feval(fun,gx,gy,gz);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);
for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem1.bas, ':', gx, gy, gz, i, MatInd, 1).*fem1.t_e_orit(:,i);
    Ibasy{i} = feEvalBas1(fem1.bas, ':', gx, gy, gz, i, MatInd, 2).*fem1.t_e_orit(:,i);
    Ibasz{i} = feEvalBas1(fem1.bas, ':', gx, gy, gz, i, MatInd, 3).*fem1.t_e_orit(:,i);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem2.bas, ':', gx, gy, gz, j, MatInd, 1).*fem2.t_e_orit(:,j);
    Jbasy{j} = feEvalBas2(fem2.bas, ':', gx, gy, gz, j, MatInd, 2).*fem2.t_e_orit(:,j);
    Jbasz{j} = feEvalBas2(fem2.bas, ':', gx, gy, gz, j, MatInd, 3).*fem2.t_e_orit(:,j);
end

I = reshape(repmat(fem1.g2ldof,6,1),nloc*nt,1);
J = repmat(reshape(fem2.g2ldof,dof2*nt,1),6,1);
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
S = sparse(I(ID),J(ID),X(ID),size(fem1.gdof,1),size(fem2.gdof,1));
