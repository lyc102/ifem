function M = globMatrixNedIFE3DFace(fun1,fun2,fem1,fem2,femIF,d1,j1,d2,j2)

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
IfeEvalBas1 = @evalNed1IFEBas3D;
IfeEvalBas2 = @evalNed1IFEBas3D;

%% 1. Matrix on Internal Interface Faces
A = femIF.area; nm = femIF.normal;
gw = femIF.gw; gx = femIF.gx; gy = femIF.gy; gz = femIF.gz;
tIntfLID = femIF.tIntfLID; tIntfRID = femIF.tIntfRID;
nt = size(femIF.tL,1); % NOT # of interface elements, but quadrature element
ntot = 4*nt*dof1*dof2;
I = zeros(ntot,1); J = zeros(ntot,1); X = zeros(ntot,1);
t_e_orit = fem1.t_e_orit;

coef1 = feval(fun1,gx,gy,gz);
coef2 = feval(fun2,gx,gy,gz);

IbasLx = cell(dof1,1); IbasRx = cell(dof1,1);
IbasLy = cell(dof1,1); IbasRy = cell(dof1,1);
IbasLz = cell(dof1,1); IbasRz = cell(dof1,1);
if d1 == 0
    for i = 1:dof1
        IbasLxtmp = IfeEvalBas1(femIF.basL(:,:,i), gx, gy, gz, 0, 1).*t_e_orit(tIntfLID,i);
        IbasLytmp = IfeEvalBas1(femIF.basL(:,:,i), gx, gy, gz, 0, 2).*t_e_orit(tIntfLID,i);
        IbasLztmp = IfeEvalBas1(femIF.basL(:,:,i), gx, gy, gz, 0, 3).*t_e_orit(tIntfLID,i);
        
        IbasRxtmp = IfeEvalBas1(femIF.basR(:,:,i), gx, gy, gz, 0, 1).*t_e_orit(tIntfRID,i);
        IbasRytmp = IfeEvalBas1(femIF.basR(:,:,i), gx, gy, gz, 0, 2).*t_e_orit(tIntfRID,i);
        IbasRztmp = IfeEvalBas1(femIF.basR(:,:,i), gx, gy, gz, 0, 3).*t_e_orit(tIntfRID,i);
        
        IbasLx{i} = IbasLytmp.*nm(:,3) - IbasLztmp.*nm(:,2);
        IbasLy{i} = IbasLztmp.*nm(:,1) - IbasLxtmp.*nm(:,3);
        IbasLz{i} = IbasLxtmp.*nm(:,2) - IbasLytmp.*nm(:,1);
        
        IbasRx{i} = IbasRytmp.*nm(:,3) - IbasRztmp.*nm(:,2);
        IbasRy{i} = IbasRztmp.*nm(:,1) - IbasRxtmp.*nm(:,3);
        IbasRz{i} = IbasRxtmp.*nm(:,2) - IbasRytmp.*nm(:,1);
        
    end
elseif d1 == 1
    for i = 1:dof1
        IbasLx{i} = IfeEvalBas1(femIF.basL(:,:,i), gx, gy, gz, 1, 1).*t_e_orit(tIntfLID,i);
        IbasLy{i} = IfeEvalBas1(femIF.basL(:,:,i), gx, gy, gz, 1, 2).*t_e_orit(tIntfLID,i);
        IbasLz{i} = IfeEvalBas1(femIF.basL(:,:,i), gx, gy, gz, 1, 3).*t_e_orit(tIntfLID,i);
        
        IbasRx{i} = IfeEvalBas1(femIF.basR(:,:,i), gx, gy, gz, 1, 1).*t_e_orit(tIntfRID,i);
        IbasRy{i} = IfeEvalBas1(femIF.basR(:,:,i), gx, gy, gz, 1, 2).*t_e_orit(tIntfRID,i);
        IbasRz{i} = IfeEvalBas1(femIF.basR(:,:,i), gx, gy, gz, 1, 3).*t_e_orit(tIntfRID,i);
        
    end
end

JbasLx = cell(dof2,1);  JbasRx = cell(dof2,1);
JbasLy = cell(dof2,1);  JbasRy = cell(dof2,1);
JbasLz = cell(dof2,1);  JbasRz = cell(dof2,1);
if d2 == 0
    for j = 1:dof2
        JbasLxtmp = IfeEvalBas2(femIF.basL(:,:,j), gx, gy, gz, 0, 1).*t_e_orit(tIntfLID,j);
        JbasLytmp = IfeEvalBas2(femIF.basL(:,:,j), gx, gy, gz, 0, 2).*t_e_orit(tIntfLID,j);
        JbasLztmp = IfeEvalBas2(femIF.basL(:,:,j), gx, gy, gz, 0, 3).*t_e_orit(tIntfLID,j);
        
        JbasRxtmp = IfeEvalBas2(femIF.basR(:,:,j), gx, gy, gz, 0, 1).*t_e_orit(tIntfRID,j);
        JbasRytmp = IfeEvalBas2(femIF.basR(:,:,j), gx, gy, gz, 0, 2).*t_e_orit(tIntfRID,j);
        JbasRztmp = IfeEvalBas2(femIF.basR(:,:,j), gx, gy, gz, 0, 3).*t_e_orit(tIntfRID,j);
        
        JbasLx{j} = JbasLytmp.*nm(:,3) - JbasLztmp.*nm(:,2);
        JbasLy{j} = JbasLztmp.*nm(:,1) - JbasLxtmp.*nm(:,3);
        JbasLz{j} = JbasLxtmp.*nm(:,2) - JbasLytmp.*nm(:,1);
       
        JbasRx{j} = JbasRytmp.*nm(:,3) - JbasRztmp.*nm(:,2);
        JbasRy{j} = JbasRztmp.*nm(:,1) - JbasRxtmp.*nm(:,3);
        JbasRz{j} = JbasRxtmp.*nm(:,2) - JbasRytmp.*nm(:,1);
        
    end
elseif d2 == 1
    for j = 1:dof2
        
        JbasLx{j} = IfeEvalBas2(femIF.basL(:,:,j), gx, gy, gz, 1, 1).*t_e_orit(tIntfLID,j);
        JbasLy{j} = IfeEvalBas2(femIF.basL(:,:,j), gx, gy, gz, 1, 2).*t_e_orit(tIntfLID,j);
        JbasLz{j} = IfeEvalBas2(femIF.basL(:,:,j), gx, gy, gz, 1, 3).*t_e_orit(tIntfLID,j);
        
        JbasRx{j} = IfeEvalBas2(femIF.basR(:,:,j), gx, gy, gz, 1, 1).*t_e_orit(tIntfRID,j);
        JbasRy{j} = IfeEvalBas2(femIF.basR(:,:,j), gx, gy, gz, 1, 2).*t_e_orit(tIntfRID,j);
        JbasRz{j} = IfeEvalBas2(femIF.basR(:,:,j), gx, gy, gz, 1, 3).*t_e_orit(tIntfRID,j);
    end
end

% Jump and Average
if j1 == 1
    for i = 1:dof1
        IbasLx{i} = IbasLx{i}; IbasRx{i} = -1*IbasRx{i};
        IbasLy{i} = IbasLy{i}; IbasRy{i} = -1*IbasRy{i};
        IbasLz{i} = IbasLz{i}; IbasRz{i} = -1*IbasRz{i};
    end
elseif j1 == 0
    for i = 1:dof1
        IbasLx{i} = 1/2*IbasLx{i}; IbasRx{i} = 1/2*IbasRx{i};
        IbasLy{i} = 1/2*IbasLy{i}; IbasRy{i} = 1/2*IbasRy{i};
        IbasLz{i} = 1/2*IbasLz{i}; IbasRz{i} = 1/2*IbasRz{i};
    end
end
if j2 == 1
    for j = 1:dof2
        JbasLx{j} = JbasLx{j};  JbasRx{j} = -1*JbasRx{j};
        JbasLy{j} = JbasLy{j};  JbasRy{j} = -1*JbasRy{j};
        JbasLz{j} = JbasLz{j};  JbasRz{j} = -1*JbasRz{j};
    end
elseif j2 == 0
    for j = 1:dof2
        JbasLx{j} = 1/2*JbasLx{j}; JbasRx{j} = 1/2*JbasRx{j};
        JbasLy{j} = 1/2*JbasLy{j}; JbasRy{j} = 1/2*JbasRy{j};
        JbasLz{j} = 1/2*JbasLz{j}; JbasRz{j} = 1/2*JbasRz{j};
    end
end

ind = 0;
for i = 1:dof1
    for j = 1:dof2
        I(ind+1:ind+4*nt) = [femIF.tL(:,i);femIF.tL(:,i);femIF.tR(:,i);femIF.tR(:,i)];
        J(ind+1:ind+4*nt) = [femIF.tL(:,j);femIF.tR(:,j);femIF.tL(:,j);femIF.tR(:,j)];
        
        X1 = A.*sum((((coef1.*IbasLx{i}).*(coef2.*JbasLx{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasLy{i}).*(coef2.*JbasLy{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasLz{i}).*(coef2.*JbasLz{j})).*gw'),2);
        X2 = A.*sum((((coef1.*IbasLx{i}).*(coef2.*JbasRx{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasLy{i}).*(coef2.*JbasRy{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasLz{i}).*(coef2.*JbasRz{j})).*gw'),2);
        X3 = A.*sum((((coef1.*IbasRx{i}).*(coef2.*JbasLx{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasRy{i}).*(coef2.*JbasLy{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasRz{i}).*(coef2.*JbasLz{j})).*gw'),2);
        X4 = A.*sum((((coef1.*IbasRx{i}).*(coef2.*JbasRx{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasRy{i}).*(coef2.*JbasRy{j})).*gw'),2) +...
             A.*sum((((coef1.*IbasRz{i}).*(coef2.*JbasRz{j})).*gw'),2);
        X(ind+1:ind+4*nt) = [...
            X1; ...
            X2; ...
            X3; ...
            X4];
        ind = ind + 4*nt;
    end
end
M = sparse(I,J,X,size(fem1.gdof,1),size(fem2.gdof,1));

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