function M = globMatrixIFE3DFace(fun1,fun2,fem1,fem2,femIF,d1,j1,d2,j2)

%% USAGE: generate global matrix on Interface Faces    
%             e.g.  \int_{Omega^I} [coef1 Du*n] {coef2 v} dx 
%
% INPUTS:
% fun1 --- coefficient function in front of test function
% fun2 --- coefficient function in front of trial function
% fem1 --- global DoF for test function space
% fem2 --- global DoF for trial function space
% femIF --- quadrature info on interface faces
% d1 --- derivative info for test function
% d2 --- derivative info for trial function
%            d = 0: function value
%            d = 1: normal gradient value Du*n
% j1 --- jump info for test function
% j2 --- jump info for trial function
%            j = 0: average e.g. {u} = 1/2*(uL+uR)
%            j = 1: average e.g. [u] = uL-uR
%            Note: on Dirichlet boundary: [u] = {u} = u
% OUTPUTS:
% [I J X] --- triplets of the sparse matrix. 

% Last Modified: 08/07/2020 by Xu Zhang

%% 0. Initialization
dof1 = fem1.ldof; dof2 = fem2.ldof;
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

%% 1. Matrix on Internal Interface Faces
A = femIF.area; nm = femIF.normal;
gw = femIF.gw; gx = femIF.gx; gy = femIF.gy; gz = femIF.gz;
nt = size(femIF.tL,1); % NOT # of interface elements, but quadrature element
ntot = 4*nt*dof1*dof2;
I = zeros(ntot,1); J = zeros(ntot,1); X = zeros(ntot,1);

coef1 = feval(fun1,gx,gy,gz);
coef2 = feval(fun2,gx,gy,gz);

IbasL = cell(dof1,1); IbasR = cell(dof1,1);
if d1 == 0
    for i = 1:dof1
        IbasL{i} = feEvalBas1(femIF.basL(:,:,i), gx, gy, gz, [0,0,0]);
        IbasR{i} = feEvalBas1(femIF.basR(:,:,i), gx, gy, gz, [0,0,0]);
    end
elseif d1 == 1
    for i = 1:dof1
        IbasLx = feEvalBas1(femIF.basL(:,:,i), gx, gy, gz, [1,0,0]);
        IbasLy = feEvalBas1(femIF.basL(:,:,i), gx, gy, gz, [0,1,0]);
        IbasLz = feEvalBas1(femIF.basL(:,:,i), gx, gy, gz, [0,0,1]);
        
        IbasRx = feEvalBas1(femIF.basR(:,:,i), gx, gy, gz, [1,0,0]);
        IbasRy = feEvalBas1(femIF.basR(:,:,i), gx, gy, gz, [0,1,0]);
        IbasRz = feEvalBas1(femIF.basR(:,:,i), gx, gy, gz, [0,0,1]);
        
        IbasL{i} = IbasLx.*nm(:,1) + IbasLy.*nm(:,2) + IbasLz.*nm(:,3);
        IbasR{i} = IbasRx.*nm(:,1) + IbasRy.*nm(:,2) + IbasRz.*nm(:,3);
    end
end

JbasL = cell(dof2,1);  JbasR = cell(dof2,1);
if d2 == 0
    for j = 1:dof2
        JbasL{j} = feEvalBas2(femIF.basL(:,:,j), gx, gy, gz, [0,0,0]);
        JbasR{j} = feEvalBas2(femIF.basR(:,:,j), gx, gy, gz, [0,0,0]);
    end
elseif d2 == 1
    for j = 1:dof2
        JbasLx = feEvalBas2(femIF.basL(:,:,j), gx, gy, gz, [1,0,0]);
        JbasLy = feEvalBas2(femIF.basL(:,:,j), gx, gy, gz, [0,1,0]);
        JbasLz = feEvalBas2(femIF.basL(:,:,j), gx, gy, gz, [0,0,1]);
        JbasRx = feEvalBas2(femIF.basR(:,:,j), gx, gy, gz, [1,0,0]);
        JbasRy = feEvalBas2(femIF.basR(:,:,j), gx, gy, gz, [0,1,0]);
        JbasRz = feEvalBas2(femIF.basR(:,:,j), gx, gy, gz, [0,0,1]);
        JbasL{j} = JbasLx.*nm(:,1) + JbasLy.*nm(:,2) + JbasLz.*nm(:,3);
        JbasR{j} = JbasRx.*nm(:,1) + JbasRy.*nm(:,2) + JbasRz.*nm(:,3);
    end
end

% Jump and Average
if j1 == 1
    for i = 1:dof1
        IbasL{i} = IbasL{i}; IbasR{i} = -1*IbasR{i};
    end
elseif j1 == 0
    for i = 1:dof1
        IbasL{i} = 1/2*IbasL{i}; IbasR{i} = 1/2*IbasR{i};
    end
end
if j2 == 1
    for j = 1:dof2
        JbasL{j} = JbasL{j};  JbasR{j} = -1*JbasR{j};
    end
elseif j2 == 0
    for j = 1:dof2
        JbasL{j} = 1/2*JbasL{j}; JbasR{j} = 1/2*JbasR{j};
    end
end

ind = 0;
for i = 1:dof1
    for j = 1:dof2
        I(ind+1:ind+4*nt) = [femIF.tL(:,i);femIF.tL(:,i);femIF.tR(:,i);femIF.tR(:,i)];
        J(ind+1:ind+4*nt) = [femIF.tL(:,j);femIF.tR(:,j);femIF.tL(:,j);femIF.tR(:,j)];
        X(ind+1:ind+4*nt) = [...
            A.*sum((((coef1.*IbasL{i}).*(coef2.*JbasL{j})).*gw'),2); ...
            A.*sum((((coef1.*IbasL{i}).*(coef2.*JbasR{j})).*gw'),2); ...
            A.*sum((((coef1.*IbasR{i}).*(coef2.*JbasL{j})).*gw'),2); ...
            A.*sum((((coef1.*IbasR{i}).*(coef2.*JbasR{j})).*gw'),2)];
        ind = ind + 4*nt;
    end
end
M = sparse(I,J,X,size(fem1.p,1),size(fem2.p,1));

%% 2. Matrix on Boundary Interface Faces
if size(femIF.tB,1) ~= 0
    nmB = femIF.normalB;
    AB = femIF.areaB; gxB = femIF.gxB; gyB = femIF.gyB; gzB = femIF.gzB;
    ntB = size(femIF.tB,1); % not number of interface element, but quadrature element
    ntotB = ntB*dof1*dof2;
    IB = zeros(ntotB,1); JB = zeros(ntotB,1); XB = zeros(ntotB,1);
    
    coef1B = feval(fun1,gxB,gyB,gzB);
    coef2B = feval(fun2,gxB,gyB,gzB);
    IbasB = cell(dof1,1);
    if d1 == 0
        for i = 1:dof1
            IbasB{i} = feEvalBas1(femIF.basB(:,:,i), gxB, gyB, gzB, [0,0,0]);
        end
    elseif d1 == 1
        for i = 1:dof1
            IbasBx = feEvalBas1(femIF.basB(:,:,i), gxB, gyB, gzB, [1,0,0]);
            IbasBy = feEvalBas1(femIF.basB(:,:,i), gxB, gyB, gzB, [0,1,0]);
            IbasBz = feEvalBas1(femIF.basB(:,:,i), gxB, gyB, gzB, [0,0,1]);
            IbasB{i} = IbasBx.*nmB(:,1) + IbasBy.*nmB(:,2) + IbasBz.*nmB(:,3);
        end
    end
    
    JbasB = cell(dof2,1);
    if d2 == 0
        for j = 1:dof2
            JbasB{j} = feEvalBas2(femIF.basB(:,:,j), gxB, gyB, gzB, [0,0,0]);
        end
    elseif d2 == 1
        for j = 1:dof2
            JbasBx = feEvalBas2(femIF.basB(:,:,j), gxB, gyB, gzB, [1,0,0]);
            JbasBy = feEvalBas2(femIF.basB(:,:,j), gxB, gyB, gzB, [0,1,0]);
            JbasBz = feEvalBas2(femIF.basB(:,:,j), gxB, gyB, gzB, [0,0,1]);
            JbasB{j} = JbasBx.*nmB(:,1) + JbasBy.*nmB(:,2) + JbasBz.*nmB(:,3);
        end
    end
    
    ind = 0;
    for i = 1:dof1
        for j = 1:dof2
            IB(ind+1:ind+ntB) = femIF.tB(:,i);
            JB(ind+1:ind+ntB) = femIF.tB(:,j);
            XB(ind+1:ind+ntB) = AB.*sum((((coef1B.*IbasB{i}).*(coef2B.*JbasB{j})).*gw'),2);
            ind = ind + ntB;
        end
    end
    MB = sparse(IB,JB,XB,size(fem1.p,1),size(fem2.p,1));
    M = M + MB;
end
