%% Check rate of convergence for 3D Hcurl interface problem
%     curl(A curl u) + B u  = f,    x\in \Omega
%      where A and B are piecewise constants on Omega^+ and Omega^-.
%
% Domain: Rectangular domain: [xmin,xmax] X [ymin,ymax] X [zmin,zmax]
% Mesh: Cartesian triangular mesh.
% Method: FE
% for this method, need to mannualy import a mesh data

%% Geometry and Boundary Conditions
%clear
%close all
%clc

%path(pathdef)
addpath(genpath(pwd),'-begin');
rmpath(genpath('./.git'));
rmpath(genpath('./docs'));
savepath;

domain = [-1,1,-1,1,-1,1];
bc = [1,1,1,1,1,1]; % Dirichelet BC

%% Finite Element Type
femtype = '1st kind Nedelec';
disp(['FEM Type =  Conforming ', femtype]);

%% Initial Partition
nx0 = 10;
ny0 = nx0;
nz0 = nx0;

%% Task
showErr = 0;
showMesh = 0;
computErr = 1;

%% PDE
test = 3;
switch test
    case 0
        pde = poissonNonPoly3D;
    case 1 % circular interface
        r = pi/5; bm = 1; bp = 1; am = 1; ap = 1;
        x0 = 0; y0 = 0; z0 = 0; a11 = 1; a12 = 1; a = 1;
        pde = elli3DcircIntf2(am,ap,bm,bp,r,x0,y0,z0,a11,a12,a);
    case 2
        r = pi/5; bm = 1; bp = 1;
        x0 = 0; y0 = 0; z0 = 0; coef = 1;
        pde = constant_fun(bm,bp,r,x0,y0,z0,coef);
    case 3
        r = pi/5; bm = 1; bp = 1; am = 1; ap = 1;
        x0 = 0; y0 = 0; z0 = 0; coef = 1;
        pde = linear_fun(am,ap,bm,bp,r,x0,y0,z0,coef);
end

%% Max Iteration

%% 1. Generate Mesh
disp('*******************************************************************************')
disp('Import mesh');
load('BoxTorusMesh.mat')
mesh = enrichMesh3D(mesh,0); % Mesh detail level = 1 (for IFE).

%% 2. Generate FEM DoF
fem = genNedFEM3D(mesh,bc);
disp(['number of DoF =  ', int2str(length(fem.p))]);

%% 3. Assemble Matrix
disp(' '); disp('Start Assembling Matrix');
S = globMatrixNedFit3D(am,ap,1,mesh,fem,fem);
M = globMatrixNedFit3D(bm,bp,0,mesh,fem,fem);
rhsF1 = globNedFitRHS3D(pde.fm1,pde.fp1, mesh, fem, 0, 1);
rhsF2 = globNedFitRHS3D(pde.fm2,pde.fp2, mesh, fem, 0, 2);
rhsF3 = globNedFitRHS3D(pde.fm3,pde.fp3, mesh, fem, 0, 3);
rhsF = rhsF1 + rhsF2 + rhsF3;
Atotal = S + M;
tu1 = sum(feval(pde.exactu1,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tu2 = sum(feval(pde.exactu2,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tu3 = sum(feval(pde.exactu3,fem.gex,fem.gey,fem.gez).*fem.gew,2);
tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
tgt = tgt./sum(tgt.^2,2).^(1/2);
tu = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);

Ndof = size(mesh.e,1);
bdidx = zeros(Ndof,1);
isBdEdge = true(Ndof,1);
isBdEdge(fem.mapper) = false;
bdidx(isBdEdge) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
A = T*Atotal*T + Tbd;
ub = tu;
ub(fem.mapper) = 0;
rhsB = Atotal*ub;
f = rhsF - rhsB;
f(isBdEdge) = tu(isBdEdge);

%% 3. Solve the linear system Au = f
option.outsolver = 'cg';
tID1 = (mesh.tLoc == 1);
tID2 = (mesh.tLoc == 2);
e1tmp = unique(reshape(mesh.t_e(tID1,:),[],1));
e2tmp = unique(reshape(mesh.t_e(tID2,:),[],1));
eInt = intersect(e1tmp,e2tmp);
%eIntp1 = mesh.p(mesh.e(eInt,1),:);
%eIntp2 = mesh.p(mesh.e(eInt,2),:);
eid1 = setdiff(e1tmp,eInt);
eid2 = setdiff(e2tmp,eInt);
alpha = am*ones(size(mesh.e,1),1);
alpha(eid2) = ap;
alpha(eInt) = (am+ap)/2;
beta = bm*ones(size(mesh.e,1),1);
beta(eid2) = bp;
beta(eInt) = (bm+bp)/2;
option.alpha = alpha;
option.beta = beta;
option.solver = 'amg';
edge = mesh.e;
option.isBdEdge = isBdEdge;
option.smoother = 'BD';
NEdof = size(A,1);
option.blklevel = 0;
option.blkId = NEdof;
option.fact = 'chol';
[x,info] = amgMaxwellinterface(A,f,mesh.p,edge,option);
uh = x;

%% 4. Postprocess: Calculating Errors
tic
errND = max(abs(uh - tu)); % Error on nodes
err.nd = errND; err.inf = 0; err.l2 = 0; err.h1 = 0;
if computErr == 1
    disp(' ');  disp('Start computing error in L2 norm');
    eNorm = 'L2'; disp(['Start computing error in ',eNorm,'  norm']);
    [errL2,errL2K1,errL2K2,errL2K3] = getCurlErr3D(uh, pde, fem, eNorm);
    
    eNorm = 'Curl'; disp(['Start computing error in ',eNorm,' norm']);
    [errCurl,errCurl1,errCurl2,errCurl3] = getCurlErr3D(uh, pde, fem, eNorm);
    err.l2 = errL2;
    err.curl = errCurl;
end
time(6) = toc;
time(7) = 1e6*sum(time(2:6))/length(mesh.t);

%% 5: Output

disp(' ')
disp('Errors')
disp('Node   L2 norm     H1 norm')
formatSpec = '%6.4e %6.4e  %6.4e\n';
fprintf(formatSpec, err.nd, err.l2, err.curl)

err0 = err; h0 = h;

disp(' '); disp('CPU Time')
disp('   N     Mesh     FEM      Matrix   Solve    Error    Time/1M cell')
formatSpec = '%4i  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f   %7.2f\n';
fprintf(formatSpec, time)

%% 6. Plot Solution and Error
if showErr == 1
end