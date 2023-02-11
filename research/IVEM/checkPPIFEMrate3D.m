%% Check rate of convergence for 3D H1 interface problem
%     -div(A grad u)  = f,    x\in \Omega
%      where A is a piecewise constant on Omega^+ and Omega^-.
%
% Domain: Rectangular domain: [xmin,xmax] X [ymin,ymax] X [zmin,zmax]
% Mesh: Cartesian triangular mesh.
% Method: PP-IFE


%% Geometry and Boundary Conditions
% clear
% close all
%clc

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
test = 6;
switch test
    case 0
        pde = poissonNonPoly3D;
    case 1 % circular interface
        r = pi/4; bm = 1; bp = pi/(2*r^2);
        %r = sqrt(pi/2); bm = 1; bp = pi/(2*r^2);domain = [-2,2,-2,2,-2,2];
        x0 = 0; y0 = 0; z0 = 0; rx = r; ry = r; rz = r;
        pde = elli3DcircIntf(bm,bp,r,x0,y0,z0,rx,ry,rz);
    case 2 % orthotorus interface
        domain = [-1.2,1.2,-1.2,1.2,-1.2,1.2];
        bm = 1; bp = 1;
        rx = 1; ry = 0.075; rz = 3;
        pde = elli3DorthocircIntf(bm,bp,rx,ry,rz);
    case 3 % line interface
        bm = 1; bp = 1;
        rx = 1; ry = 10^(-4); rz = -1; cx = 10^(-4); cy = 0; cz = 0;
        cm = 1; cp = 1; a = 1;
        pde = elli3DlinIntf(bm,bp,cx,cy,cz,rx,ry,rz,a,cm,cp);
    case 4 % line interface but with zero boundary conditions
        bm = 1; bp = 10; delta = pi/100; % delta can not be -1,1
        pde = elli3DlinIntf2(bm,bp,delta);
    case 5 % circular interface
        r = pi/4; bm = 1; bp = 10;domain = [-1,1,-1,1,-1,1];
        x0 = 0; y0 = 0; z0 = 0;
        pde = elli3DcircIntf5(bm,bp,r,x0,y0,z0);
    case 6 % circular interface
        bm = 1; bp = 1;domain = [-1.3,1.3,-1.3,1.3,-1.3,1.3];
        x1 = -0.3; y1 = 0; z1 = 0; r11 = pi/5; r12 = 0.2;
        x2 =  0.3; y2 = 0; z2 = 0; r21 = pi/5; r22 = 0.2;
        pde = elli3DtorusTwin(bm,bp,x1,y1,z1,r11,r12,x2,y2,z2,r21,r22);
end
%% PPIFEM Type
PPtype = 'S';
disp(['PPIFEM Type = ', PPtype]);
sig = max(bm,bp)/min(bm,bp);

%% Max Iteration
maxIt = 10;
time = zeros(9,maxIt);

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
    mesh = enrichMesh3D(mesh,2); % Mesh detail level = 2 (for PPIFE).
    disp(['number of element =  ', int2str(length(mesh.t))]);
    mesh = genIntfMesh3D(mesh,pde.intf);
    disp(['number of interface element =  ', int2str(-min(mesh.tLoc)),...
        ', is ', num2str(100*-min(mesh.tLoc)/length(mesh.t)), '% of all elements']);
    if showMesh == 1
        tetramesh(mesh.t,mesh.p,ones(size(mesh.t)));
    end
    time(2,i) = toc;
    
    %% 2. Generate FEM DoF
    tic
    fem = genFEM3D(mesh,femtype,bc);
    disp(['number of DoF =  ', int2str(length(fem.p))]);
    time(3,i) = toc;
    tic
    femI = genP1IFEM3D(mesh,fem,bm,bp);
    time(4,i) = toc;
    tic
    femIF = genP1IFEM3DFace(mesh,fem,femI);
    time(5,i) = toc;
    
    %% 3. Assemble Matrix
    tic
    disp(' '); disp('Start Assembling Matrix');
    matrix = genMatEllPPIFE3D(pde,mesh,fem,femI,femIF,PPtype,sig);
    time(6,i) = toc;
    
    %% 4. Solve the linear system Au = f
    tic
    disp(' ')
    if strcmp(PPtype,'S')
        disp('Start Solving Linear System: using PCG with ichol precond');
        L = ichol(matrix.A);
        u = pcg(matrix.A,matrix.f,1e-8,300,L,L');
    elseif strcmp(PPtype,'N') || strcmp(PPtype,'I')
        disp('Start Solving Linear System: using GMRES with ilu precond');
        [L,U] = ilu(matrix.A,struct('milu','row'));
        u = gmres(matrix.A,matrix.f,[],1e-8,300,L,U);
    end
    time(7,i) = toc;
    tu = pde.exactu(fem.p(:,1),fem.p(:,2),fem.p(:,3));
    uh = tu; uh(fem.mapper) = u;
    
    %% 5. Postprocess: Calculating Errors
    tic
    errND = max(abs(uh - tu)); % Error on nodes
    err.nd = errND; err.inf = 0; err.l2 = 0; err.h1 = 0;
    if computErr == 1
        disp(' ');  disp('Start computing error in Inf norm');
        [errInfN,errInfI] = getErrInfIFE3D(uh,pde.exactu,mesh,femI,fem,[0,0,0]);        
        eNorm = 'L2'; disp(['Start computing error in ',eNorm,'  norm']);
        [errL2,errL2K] = getErrIFE3D(uh, pde, mesh, fem, femI,eNorm);
        eNorm = 'H1x'; disp(['Start computing error in ',eNorm,' norm']);
        [errH1x,errH1xK] = getErrIFE3D(uh, pde, mesh, fem, femI, eNorm);
        eNorm = 'H1y'; disp(['Start computing error in ',eNorm,' norm']);
        [errH1y,errH1yK] = getErrIFE3D(uh, pde, mesh, fem, femI, eNorm);
        eNorm = 'H1z'; disp(['Start computing error in ',eNorm,' norm']);
        [errH1z,errH1zK] = getErrIFE3D(uh, pde, mesh, fem, femI, eNorm);
        
        disp(' ')
        disp(['Max error on Regular Elements   = ', num2str(errInfN)])
        disp(['Max error on Interface Elements = ', num2str(errInfI)])
        err.inf = max([errND,errInfN,errInfI]);
        err.l2 = errL2;
        err.h1 = sqrt(errH1x^2+errH1y^2+errH1z^2);
    end
    time(8,i) = toc;
    time(9,i) = 1e6*sum(time(2:8,i))/length(mesh.t);
    
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
    end
    err0 = err; h0 = h;
    
    disp(' '); disp('CPU Time')
    disp('   N     Mesh     FEM      FEMI    FEMIF    Matrix   Solve    Error    Time/1M cell')
    formatSpec = '%4i  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f   %7.2f\n';
    fprintf(formatSpec, time(:,1:i))
    
    %% 6. Plot Solution and Error
    if showErr == 1
    end
   
end