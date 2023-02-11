%% Check rate of convergence for 3D Hcurl interface problem
%     curl(A curl u) + B u  = f,    x\in \Omega
%      where A and B are piecewise constants on Omega^+ and Omega^-.
%
% Domain: Rectangular domain: [xmin,xmax] X [ymin,ymax] X [zmin,zmax]
% Mesh: Cartesian triangular mesh.
% Method: FE without interface

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
test = 1;
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
maxIt = 7;
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
    mesh = enrichMesh3D(mesh,0); % Mesh detail level = 1 (for IFE).
    mesh = genIntfMesh3D(mesh,pde.intf);
    disp(['number of interface element =  ', int2str(-min(mesh.tLoc)),...
        ', is ', num2str(100*-min(mesh.tLoc)/length(mesh.t)), '% of all elements']);
    if showMesh == 1
        tetramesh(mesh.t,mesh.p,ones(size(mesh.t)));
    end
    time(2,i) = toc;
    
    %% 2. Generate FEM DoF
    tic
    fem = genNedFEM3D(mesh,bc);
    disp(['number of DoF =  ', int2str(length(fem.p))]);
    time(3,i) = toc;
    
    %% 3. Assemble Matrix
    tic
    disp(' '); disp('Start Assembling Matrix');
    matrix = genMatCurl3D(pde,mesh,fem);
    time(4,i) = toc;
    
    %% 3. Solve the linear system Au = f
    tic
    % L = ichol(matrix.A,struct('michol','on'));
    L = ichol(matrix.A);
    [u,flag,relres,iter,resvec] = pcg(matrix.A,matrix.f,1e-8,500,L,L');
    time(5,i) = toc;
    tu = matrix.tu;
    uh = tu; uh(fem.mapper) = u;
%     uh = matrix.tu; tu=uh;
    
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
    time(6,i) = toc;
    time(7,i) = 1e6*sum(time(2:6,i))/length(mesh.t);
    
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
    disp('   N     Mesh     FEM      Matrix   Solve    Error    Time/1M cell')
    formatSpec = '%4i  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f   %7.2f\n';
    fprintf(formatSpec, time(:,1:i))
    
    %% 6. Plot Solution and Error
    if showErr == 1
    end
end