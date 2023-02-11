function rhs = globRHSIFE3DFace(fun, funB, fem, femIF, d1)

%% USAGE: generate global matrix on a 3D mesh Face
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
% matrix --- global (mass, stiffness, ...) matrix.

% Last Modified: 08/02/2020 by Xu Zhang

%% 1. RHS on Boundary Interface Faces
dof = fem.ldof;
if strcmp(fem.type,'P1')||strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR')
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem.type,'P2')||strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end
rhs = zeros(length(fem.p),1);

if size(femIF.tB,1) ~= 0
    
    nmB = femIF.normalB;
    AB = femIF.areaB; gxB = femIF.gxB; gyB = femIF.gyB; gzB = femIF.gzB; gw = femIF.gw;
    nt = size(femIF.tB,1); % not number of interface element, but quadrature element
    ntot = nt*dof;
    IB = zeros(ntot,1); XB = zeros(ntot,1);
    
    coef = feval(fun,gxB,gyB,gzB);
    bval = feval(funB,gxB,gyB,gzB);
    IbasB = cell(dof,1);
    if d1 == 0
        for i = 1:dof
            IbasB{i} = feEvalBas(femIF.basB(:,:,i), gxB, gyB, gzB, [0,0,0]);
        end
    elseif d1 == 1
        for i = 1:dof
            IbasBx = feEvalBas(femIF.basB(:,:,i), gxB, gyB, gzB, [1,0,0]);
            IbasBy = feEvalBas(femIF.basB(:,:,i), gxB, gyB, gzB, [0,1,0]);
            IbasBz = feEvalBas(femIF.basB(:,:,i), gxB, gyB, gzB, [0,0,1]);
            IbasB{i} = IbasBx.*nmB(:,1) + IbasBy.*nmB(:,2) + IbasBz.*nmB(:,3);
        end
    end
    
    ind = 0;
    for i = 1:dof
        IB(ind+1:ind+nt) = femIF.tB(:,i);
        XB(ind+1:ind+nt) = AB.*sum((((IbasB{i}.*coef).*bval).*gw'),2);
        ind = ind + nt;
    end
    ID = find(XB~=0);
    rhs = sparse(IB(ID),1,XB(ID),length(fem.p),1);
end

