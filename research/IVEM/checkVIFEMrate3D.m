%% Check rate of convergence for 3D H1 interface problem
%     -div(A grad u)  = f,    x\in \Omega
%      where A is a piecewise constant on Omega^+ and Omega^-.
%
% Domain: Rectangular domain: [xmin,xmax] X [ymin,ymax] X [zmin,zmax]
% Mesh: Cartesian triangular mesh.
% Method: VIFE


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
femtype = 'P1';
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
test = 5;
switch test
    case 0
        pde = poissonNonPoly3D;
    case 1 % circular interface
        %r = pi/5; bm = 1; bp = pi/(2*r^2);
        r = sqrt(pi/2); bm = 1; bp = pi/(2*r^2);domain = [-2,2,-2,2,-2,2];
        x0 = 0; y0 = 0; z0 = 0; rx = r; ry = r; rz = r;
        pde = elli3DcircIntf(bm,bp,r,x0,y0,z0,rx,ry,rz);
    case 2 % orthotorus interface
        domain = [-1.2,1.2,-1.2,1.2,-1.2,1.2];
        bm = 1; bp = 1;
        rx = 1; ry = 0.075; rz = 3;
        pde = elli3DorthocircIntf(bm,bp,rx,ry,rz);
    case 3 % line interface
        bm = 1; bp = 10;
        rx = 1; ry = 0; rz = 0; cx = pi/10; cy = 0; cz = 0;
        cm = 1; cp = 1; a = 1;
        pde = elli3DlinIntf(bm,bp,cx,cy,cz,rx,ry,rz,a,cm,cp);
    case 4 % line interface but with zero boundary conditions
        bm = 1; bp = 10^(2); delta = pi/10; % delta can not be -1,1
        pde = elli3DlinIntf2(bm,bp,delta);
    case 5 % circular interface
        r = pi/4; bm = 1; bp = 10; domain = [-1,1,-1,1,-1,1];
        x0 = 0; y0 = 0; z0 = 0;
        pde = elli3DcircIntf5(bm,bp,r,x0,y0,z0);
    case 6 % circular interface
        bm = 1; bp = 10; domain = [-1.3,1.3,-1.3,1.3,-1.3,1.3];
        x1 = -0.3; y1 = 0; z1 = 0; r11 = pi/5; r12 = 0.2;
        x2 =  0.3; y2 = 0; z2 = 0; r21 = pi/5; r22 = 0.2;
        pde = elli3DtorusTwin(bm,bp,x1,y1,z1,r11,r12,x2,y2,z2,r21,r22);
end

%% Max Iteration
maxIt = 4;
time = zeros(9,maxIt);
error = zeros(maxIt,4);
ratio = zeros(maxIt,4);

for i = 1:maxIt
    
    %% 1. Generate Mesh
    tic
    nx = nx0 + 10*(i-1); h = (domain(2) - domain(1))/nx;
    ny = ny0 + 10*(i-1);
    nz = nz0 + 10*(i-1);
    time(1,i) = nx;
    disp(' ')
    disp('**************************************************************************************')
    disp(['Partition =  ',int2str(nx),' X ',int2str(ny),' X ',int2str(nz)]);
    disp(' ')
    
    mesh = genMesh3D(domain, nx, ny, nz);
    mesh = enrichMesh3D(mesh,1); % Mesh detail level = 2 (for PPIFE).
    disp(['number of element =  ', int2str(length(mesh.t))]);    
    mesh = genIntfMesh3D(mesh,pde.intf);
    disp(['number of interface element =  ', int2str(-min(mesh.tLoc)),...
        ', is ', num2str(100*-min(mesh.tLoc)/length(mesh.t)), '% of all elements']);
    time(2,i) = toc;
    
    tic
    fem = genFEM3D(mesh,femtype,bc);
    disp(['number of DoF =  ', int2str(length(fem.p))]);
    time(3,i) = toc;
    
    %% 2. Generate FEMI data including IFE functions and their matrices
    tic
    meshI = genIVmesh(mesh);
    time(4,i) = toc;
    length(unique(union(union(meshI.tface(:,1),meshI.tface(:,2)),meshI.tface(:,3))))
    
    tic
    option.mass = 0; % generate mass matrix
    option.rhs = 1;
    femI = genP1VIFEM3D(meshI,mesh,pde,option);
    time(5,i) = toc;
    
    %% 3. Assemble Matrix
    tic
    disp(' '); disp('Start Assembling Matrix');
    option.mass = 0; % generate mass matrix
    matrix = genMatVIFE3D(pde,mesh,fem,meshI,femI);
    time(6,i) = toc;
    
    %% 4. Solve the linear system Au = f
    tic
    disp(' ')

    disp('Start Solving Linear System: using PCG with ichol precond');
    %L = ichol(matrix.A);
    alpha = 10;
    L = ichol(matrix.A,struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    [u,flag,relres,iter,resvec] = pcg(matrix.A,matrix.f,1e-8,1000,L,L');
    
    time(7,i) = toc;
    tu = matrix.tu;
    uh = tu; %uh(matrix.mapper) = u;
    
    %% 5. Postprocess: Calculating Errors
    tic
    errND = max(abs(uh - tu)); % Error on nodes
    err.nd = errND; err.inf = 0; err.l2 = 0; err.h1 = 0;
    if computErr == 1
        disp(' ');  disp('Start computing error in Inf norm');
        eNorm = 'inf';
        errInfN = getErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);
        eNorm = 'L2'; disp(['Start computing error in ',eNorm,'  norm']);
        errL2 = getErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);
        eNorm = 'H1x'; disp(['Start computing error in ',eNorm,' norm']);
        errH1x = getErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);
        eNorm = 'H1y'; disp(['Start computing error in ',eNorm,' norm']);
        errH1y = getErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);
        eNorm = 'H1z'; disp(['Start computing error in ',eNorm,' norm']);
        errH1z = getErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);
        
        disp(' ')
        err.inf = max([errND,errInfN]);
        err.l2 = errL2;
        err.h1 = sqrt(errH1x^2+errH1y^2+errH1z^2);
    end
    time(8,i) = toc;
    time(9,i) = 1e6*sum(time(2:8,i))/length(mesh.t);
    error(i,1) = errND; error(i,2) = err.inf;
    error(i,3) = err.l2; error(i,4) = err.h1;
    
    %% 6: Display Output
    disp(' ')
    disp('Errors')
    disp('Node        Inf norm    L2 norm     H1 norm')
    formatSpec = '%6.4e  %6.4e  %6.4e  %6.4e\n';
    fprintf(formatSpec, err.nd, err.inf, err.l2, err.h1)
    
    if i > 1
        format short
        rNd = log(err0.nd/err.nd)./log(h0/h);
        rInf = log(err0.inf/err.inf)./log(h0/h);
        rL2 = log(err0.l2/err.l2)./log(h0/h);
        rH1 = log(err0.h1/err.h1)./log(h0/h);
        disp(' ')
        disp('Convergence Rate')
        disp('Node        Inf norm    L2 norm     H1 norm')
        formatSpec = '%6.4f      %6.4f      %6.4f      %6.4f\n';
        fprintf(formatSpec, rNd, rInf, rL2, rH1)
        ratio(i,:) = [rNd,rInf,rL2,rH1];
    end
    err0 = err; h0 = h;
    
    disp(' '); disp('CPU Time')
    disp('   N     Mesh     FEM      MeshI    FemI    Matrix   Solve    Error    Time/1M cell')
    formatSpec = '%4i  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f   %7.2f\n';
    fprintf(formatSpec, time(:,1:i))
    
    %% 6. Plot Solution and Error
    if showErr == 1
    end
   
end