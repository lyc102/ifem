function matrix = genMatCurlVIFE3D(pde,mesh,fem,meshI,femI,option)
%% Generate global matrices and load vector of FEM for 3D ellipic eq
%     -div(A grad u)  = f,    x\in \Omega
% INPUTS:
% pde --- given data function from equation, e.g. 
%         pde.A --- diffusion coefficient
%         pde.f --- right hand side function
%         pde.gD --- Dirichlet boundary value function 
%         pde.one --- constant function 1.
% mesh --- mesh structure. 
% fem --- global degree of freedom of FEM 
%
% OUTPUTS:
% matrix.S --- stiffness matrix (w/o boundary condition)
% matrix.A --- final FEM matrix (after boundary condition)
% matrix.rhsF --- load vector (w/o boundary condition)
% matrix.f --- final RHS matrix (after boundary condition)

% Last Modified: 08/07/2020 by Xu Zhang

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g2ldofNint = femI.g2ldofNint;
node = meshI.node;
gdof = femI.gdof;
NEdof = size(gdof,1);

% %% 1. Stiffness Matrix
feEvalBas1 = @EvalNed1Bas3D;
feEvalBas2 = @EvalNed1Bas3D;

dof1 = 6; dof2 = 6; nloc = dof1*dof2; 
ntID = find(mesh.tLoc > 0); ntN = length(ntID);

%% 1. Matrix on noninterface elements
MatInd = 1;
AN = fem.area(ntID); 
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); gw = fem.gw;
XN = zeros(nloc*ntN, 1);

coefN = feval(pde.A,gxN,gyN,gzN);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);
t_e_orit = fem.t_e_orit(ntID,:);

for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 1).*t_e_orit(:,i);
    Ibasy{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 2).*t_e_orit(:,i);
    Ibasz{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 3).*t_e_orit(:,i);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 1).*t_e_orit(:,j);
    Jbasy{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 2).*t_e_orit(:,j);
    Jbasz{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 3).*t_e_orit(:,j);
end

IN = reshape(repmat(g2ldofNint(:,1:6),6,1),nloc*ntN,1);
JN = repmat(reshape(g2ldofNint(:,1:6),dof2*ntN,1),6,1);
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
SN = sparse(IN(ID),JN(ID),XN(ID),NEdof,NEdof);

%%%%%%

MatInd = 0;
AN = fem.area(ntID); 
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); gw = fem.gw;
XN = zeros(nloc*ntN, 1);

coefN = feval(pde.B,gxN,gyN,gzN);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);

for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 1).*t_e_orit(:,i);
    Ibasy{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 2).*t_e_orit(:,i);
    Ibasz{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 3).*t_e_orit(:,i);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 1).*t_e_orit(:,j);
    Jbasy{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 2).*t_e_orit(:,j);
    Jbasz{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 3).*t_e_orit(:,j);
end

IN = reshape(repmat(g2ldofNint(:,1:6),6,1),nloc*ntN,1);
JN = repmat(reshape(g2ldofNint(:,1:6),dof2*ntN,1),6,1);
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
MN = sparse(IN(ID),JN(ID),XN(ID),NEdof,NEdof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dof = 6;  nloc = dof;
X = zeros(nloc*ntN, 1);

fN1 = feval(pde.f1,gxN,gyN,gzN);
fN2 = feval(pde.f2,gxN,gyN,gzN);
fN3 = feval(pde.f3,gxN,gyN,gzN);
ind = 0;
I = reshape(g2ldofNint(:,1:6),nloc*ntN,1);
for i = 1:dof
    ibas1 = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, 0, 1);
    ibas2 = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, 0, 2);
    ibas3 = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, 0, 3);
    X(ind+1:ind+ntN) = AN.*sum((ibas1.*fN1+ibas2.*fN2+ibas3.*fN3).*gw',2).*t_e_orit(:,i);
    ind = ind + ntN;
end
rhsN = sparse(I,1,X,NEdof,1);

if isfield(option,'HM')
    if option.HM == 1
        Atotal = SN - pde.kappa*MN + femI.KI - pde.kappa*femI.SI;
    end
else
    Atotal = SN + MN + femI.KI + femI.SI;
end



rhs = rhsN + femI.b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %% 3. Dirichlet Boundary Conditions

% non-interface edges
vSign = meshI.vSign;
eSign = sum(vSign(femI.gdof),2)/2;
eIdm = find(eSign==-1/2); eIdp = find(eSign==1/2);
[gew,gex,gey,gez] = gaussPedge(node(gdof(:,1),:),node(gdof(:,2),:),2);
% tu1 = sum(feval(pde.exactu1,gex,gey,gez).*gew,2);
% tu2 = sum(feval(pde.exactu2,gex,gey,gez).*gew,2);
% tu3 = sum(feval(pde.exactu3,gex,gey,gez).*gew,2);
tu1 = sum(feval(pde.exactu1,gex,gey,gez).*gew,2);
tu2 = sum(feval(pde.exactu2,gex,gey,gez).*gew,2);
tu3 = sum(feval(pde.exactu3,gex,gey,gez).*gew,2);
tu1(eIdm) = sum(feval(pde.um1,gex(eIdm,:),gey(eIdm,:),gez(eIdm,:)).*gew(eIdm,:),2);
tu1(eIdp) = sum(feval(pde.up1,gex(eIdp,:),gey(eIdp,:),gez(eIdp,:)).*gew(eIdp,:),2);
tu2(eIdm) = sum(feval(pde.um2,gex(eIdm,:),gey(eIdm,:),gez(eIdm,:)).*gew(eIdm,:),2);
tu2(eIdp) = sum(feval(pde.up2,gex(eIdp,:),gey(eIdp,:),gez(eIdp,:)).*gew(eIdp,:),2);
tu3(eIdm) = sum(feval(pde.um3,gex(eIdm,:),gey(eIdm,:),gez(eIdm,:)).*gew(eIdm,:),2);
tu3(eIdp) = sum(feval(pde.up3,gex(eIdp,:),gey(eIdp,:),gez(eIdp,:)).*gew(eIdp,:),2);
tgt = node(gdof(:,2),:) - node(gdof(:,1),:);
tgt = tgt./sum(tgt.^2,2).^(1/2);
tu = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);

bcind = fem.bcind;
[bc,mapper] = boundaryEdge3D(node,gdof,bcind);
bdidx = zeros(NEdof,1); 
isBdEdge = true(NEdof,1);
isBdEdge(mapper) = false;
bdidx(isBdEdge) = 1;
Tbd = spdiags(bdidx,0,NEdof,NEdof);
T = spdiags(1-bdidx,0,NEdof,NEdof);
A = T*Atotal*T + Tbd;

ub = tu;
ub(mapper) = 0;
rhsB = Atotal*ub;
f = rhs - rhsB;
f(isBdEdge) = tu(isBdEdge);

%% Outputs
matrix = struct('A', A, 'f', f,'tu',tu,'isBdEdge',isBdEdge);