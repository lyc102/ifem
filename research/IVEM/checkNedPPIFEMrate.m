%% Check rate of convergence for 3D Hcurl interface problem
%     curl(A curl u) + B u  = f,    x\in \Omega
%      where A and B are piecewise constants on Omega^+ and Omega^-.
%
% Domain: Rectangular domain: [xmin,xmax] X [ymin,ymax] X [zmin,zmax]
% Mesh: Cartesian triangular mesh.
% Method: PP-IFE (It does not lead to optimal convergence)

%% Geometry and Boundary Conditions
% clear
% close all
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
test = 6;
switch test
    case 0
        pde = poissonNonPoly3D;
    case 1 % circular interface
        r = pi/5; bm = 1; bp = 20; am = 1; ap = 20;
        x0 = 0; y0 = 0; z0 = 0; a11 = 1; a12 = 1; a = 1;
        pde = elli3DcircIntf2(am,ap,bm,bp,r,x0,y0,z0,a11,a12,a);
    case 2
        r1 = 1; r2 = 0; r3 = 0; x0 = pi/10; y0 = pi/10; z0 = pi/10;
        bm = 10; bp = 1; am = 10; ap = 1;
        pde = ConstFun(am,ap,bm,bp,r1,r2,r3,x0,y0,z0);
    case 3
        r1 = 2; r2 = -1; r3 = -1; x0 = pi/10; y0 = 0; z0 = 0;
        bm = 1; bp = 10; am = 1; ap = 10;
        pde = LinearFun(am,ap,bm,bp,r1,r2,r3,x0,y0,z0);
    case 5
        r1 = pi/5; r2 = pi/2; bm = 1; bp = 10; am = 1; ap = 20;
        n2 = 20; n1 = n2*(r2^2-r1^2);
        pde = elli3DcircIntf3(am,ap,bm,bp,r1,r2,n1,n2);
    case 6 % circular interface
        bm = 1; bp = 100; am = 1; ap = 200;
        domain = [-1.3,1.3,-1.3,1.3,-1.3,1.3];
        x1 = -0.3; y1 = 0; z1 = 0; r11 = pi/5; r12 = 0.2;
        x2 =  0.3; y2 = 0; z2 = 0; r21 = pi/5; r22 = 0.2;
        pde = elli3DtorusTwinHcurl(bm,bp,am,ap,x1,y1,z1,r11,r12,x2,y2,z2,r21,r22);
end

%% PPIFEM Type
PPtype = 'S';
disp(['PPIFEM Type = ', PPtype]);
sig = 10*max(bm,bp)/min(bm,bp);
PC='BPCG';

%% Max Iteration
maxIt = 8;
time = zeros(7,maxIt);

for i = 1:maxIt
    
    %% 1. Generate Mesh
    tic
    nx = nx0 + 10*(i-1); h = (domain(2) - domain(1))/nx;
    ny = ny0 + 10*(i-1);
    nz = nz0 + 10*(i-1);
    time(1,i) = nx;
    disp(' ')
    disp('*******************************************************************************')
    disp(['Partition =  ',int2str(nx),' X ',int2str(ny),' X ',int2str(nz)]);
    disp(' ')
    
    mesh = genMesh3D(domain, nx, ny, nz);
    mesh = enrichMesh3D(mesh,2); % Mesh detail level = 1 (for IFE).
    mesh = genIntfMesh3D(mesh,pde.intf);
    
        %%%%% reoder the index of edges such that edges of interace elements
    %%%%% appear after others
    Explevel = 2; % it seems not working
    levelCount = 0;
    BasicEdgeNum = size(mesh.e,1);
    % level 1
    if Explevel>0
        OldedgeID = unique(reshape(mesh.t_e(mesh.tLoc<0,:),[],1));
        BasicEdge = true(size(mesh.e,1),1); BasicEdge(OldedgeID) = false;
        BasicEdgeNum = sum(BasicEdge);
        TotalOldEdge = zeros(size(mesh.e,1),1);
        TotalOldEdge(BasicEdge) = 1:BasicEdgeNum;
        TotalOldEdge(OldedgeID) = BasicEdgeNum+1:size(mesh.e,1);
        mesh.t_e = TotalOldEdge(mesh.t_e);
        mesh.f_e = TotalOldEdge(mesh.f_e);
        [~,TotalOldEdgeSortID] = sort(TotalOldEdge);
        mesh.eLoc = mesh.eLoc(TotalOldEdgeSortID);
        mesh.e = mesh.e(TotalOldEdgeSortID,:);
        levelCount = levelCount + 1;
    end
    tId = mesh.tLoc<0;
    while (levelCount<Explevel)
        fId = unique(reshape(mesh.t_f(tId,:),[],1));
        tId = unique(reshape(mesh.f_t(fId,:),[],1));
        tId = tId(tId>0);
        OldedgeID = unique(reshape(mesh.t_e(tId,:),[],1));
        BasicEdge = true(size(mesh.e,1),1); BasicEdge(OldedgeID) = false;
        BasicEdgeNum = sum(BasicEdge);
        TotalOldEdge = zeros(size(mesh.e,1),1);
        TotalOldEdge(BasicEdge) = 1:BasicEdgeNum;
        TotalOldEdge(OldedgeID) = BasicEdgeNum+1:size(mesh.e,1);
        mesh.t_e = TotalOldEdge(mesh.t_e);
        mesh.f_e = TotalOldEdge(mesh.f_e);
        [~,TotalOldEdgeSortID] = sort(TotalOldEdge);
        mesh.eLoc = mesh.eLoc(TotalOldEdgeSortID);
        mesh.e = mesh.e(TotalOldEdgeSortID,:);
        levelCount = levelCount + 1;
    end
    %%%%%
    
    disp(['number of interface element =  ', int2str(-min(mesh.tLoc)),...
        ', is ', num2str(100*-min(mesh.tLoc)/length(mesh.t)), '% of all elements']);
    time(2,i) = toc;
    
    %% 2. Generate FEM DoF
    tic
    fem = genNedFEM3D(mesh,bc);
    disp(['number of DoF =  ', int2str(length(fem.p))]);
    time(3,i) = toc;
    tic
    femI = genNed1IFEM3D(mesh,fem,bm,bp,am,ap);
    femIF = genNedIFEM3DFace(mesh,fem,femI);
    time(4,i) = toc;
    
    %% 3. Assemble Matrix
    tic
    disp(' '); disp('Start Assembling Matrix');
    matrix = genMatCurlPPIFE3D(pde,mesh,fem,femI,femIF,PPtype,sig);
    time(5,i) = toc;
    
    %% 3. Solve the linear system Au = f
    tic
    if strcmp(PC,'BPCG')
        tic
        disp('Start Solving Linear System: using PCG with ichol precond');
        SolverMaxIter = 10000;
        x0 = zeros(size(matrix.f));
        D = diag(matrix.A(1:BasicEdgeNum,1:BasicEdgeNum));
        M = matrix.A(BasicEdgeNum+1:end,BasicEdgeNum+1:end);
        [L,U,P,Q] = lu(M);
        [x,info] = blockpcg(matrix.A,matrix.f,x0,1e-8,SolverMaxIter,D,L,U,P,Q);
        disp(['iteration number of PCG =',int2str(info.itStep)]);
        time(6,i) = toc;
        uh = x;
        time(6,i) = toc;
    elseif strcmp(PC,'Bmg')
        option.outsolver = 'cg';
        tID1 = (mesh.tLoc == 1);
        tID2 = (mesh.tLoc == 2);
        e1tmp = unique(reshape(mesh.t_e(tID1,:),[],1));
        e2tmp = unique(reshape(mesh.t_e(tID2,:),[],1));
        eInt = intersect(e1tmp,e2tmp);
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
        option.isBdEdge = matrix.isBdEdge;
        option.smoother = 'BD';
        NEdof = size(matrix.A,1);
        option.blklevel = 0;
        option.blkId = NEdof;
        option.fact = 'chol';
        [x,info] = amgMaxwellinterface(matrix.A,matrix.f,mesh.p,edge,option);
        uh = x;
    end
    
    time(6,i) = toc;
    tu = matrix.tu;
%     uh = matrix.tu; tu=uh;
    
    %% 4. Postprocess: Calculating Errors
    tic
    errND = max(abs(uh - tu)); % Error on nodes
    err.nd = errND; err.inf = 0; err.l2 = 0; err.h1 = 0;
    if computErr == 1
        disp(' ');  disp('Start computing error in L2 norm');
        eNorm = 'L2'; disp(['Start computing error in ',eNorm,'  norm']);
        [errL2,errL2K1,errL2K2,errL2K3] = getCurlErrIFE3D(uh, pde, mesh, fem, femI, eNorm); 
        
        eNorm = 'Curl'; disp(['Start computing error in ',eNorm,' norm']);
        [errCurl,errCurl1,errCurl2,errCurl3] = getCurlErrIFE3D(uh, pde, mesh, fem, femI, eNorm);
        err.l2 = errL2;
        err.curl = errCurl;
    end
    time(7,i) = toc;
    time(8,i) = 1e6*sum(time(2:6,i))/length(mesh.t);
    
    %% 5: Output
    
    disp(' ')
    disp('Errors')
    disp('Node   L2 norm     H1 norm')
    formatSpec = '%6.4e %6.4e  %6.4e\n';
    fprintf(formatSpec, err.nd, err.l2, err.curl)
    
    if i > 1
        format short
        rNd = log(err0.nd/err.nd)./log(h0/h);
        rL2 = log(err0.l2/err.l2)./log(h0/h);
        rcurl = log(err0.curl/err.curl)./log(h0/h);
        disp(' ')
        disp('Convergence Rate')
        disp('Node     L2 norm     Curl norm')
        formatSpec = '%6.4f    %6.4f      %6.4f\n';
        fprintf(formatSpec, rNd, rL2, rcurl)
    end
    err0 = err; h0 = h;
    
    disp(' '); disp('CPU Time')
    disp('   N     Mesh     FEM   FEMI   Matrix   Solve    Error    Time/1M cell')
    formatSpec = '%4i  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f   %7.2f\n';
    fprintf(formatSpec, time(:,1:i))
    
    %% 6. Plot Solution and Error
    if showErr == 1
    end
end