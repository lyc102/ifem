%% Check rate of convergence for 3D Hcurl interface problem
%     curl(A curl u) + B u  = f,    x\in \Omega
%      where A and B are piecewise constants on Omega^+ and Omega^-.
%
% Domain: Rectangular domain: [xmin,xmax] X [ymin,ymax] X [zmin,zmax]
% Mesh: Cartesian triangular mesh.
% Method: VIFE

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
%NN = [11,21,31,41,51,61,71];
NN = [10,20,30,40,50,60,70];
%NN = [61,70];

%% Task
showErr = 0;
computErr = 1;

%% PDE
PDEoption=[];
test = 5;
switch test
    case 0
        pde = poissonNonPoly3D;
    case 1 % circular interface
        r = pi/5; bm = 1; bp = 10; am = 1; ap = 10;
        x0 = 0; y0 = 0; z0 = 0; a11 = 1; a12 = 1; a = 1;
        pde = elli3DcircIntf2(am,ap,bm,bp,r,x0,y0,z0,a11,a12,a);
    case 2
        r1 = 1; r2 = -1; r3 = -1; x0 = pi/5; y0 = pi/3; z0 = pi/4;
        bm = 100; bp = 1; am = 1; ap = 1;
        pde = ConstFun(am,ap,bm,bp,r1,r2,r3,x0,y0,z0);
    case 3 % used to test the robustness w.r.t small subelements
        r1 = 1; r2 = 0; r3 = 0; x0 = 0.05/10^4; y0 = 0; z0 = 0;
        bm = 1; bp = 10; am = 1; ap = 10;
        pde = LinearFun(am,ap,bm,bp,r1,r2,r3,x0,y0,z0);
    case 4
        x0 = 0.1;
        bm = 1; bp = 10; am = 1; ap = 10; % am = bm, ap = bp
        pde = SinFlat(am,ap,bm,bp,x0);
    case 5  % ZJ's example
        r1 = pi/4; r2 = pi/2; 
        bm = 1; bp = 10; am = 1; ap = 10;
        n2 = 20; n1 = n2*(r2^2-r1^2);
        pde = elli3DcircIntf3(am,ap,bm,bp,r1,r2,n1,n2);
    case 6
        r = pi/5; x0 = 0; y0 = 0; z0 = 0;
        bm = 1; bp = 10; am = 1; ap = 10;
        pde = LinearFun2(am,ap,bm,bp,r,x0,y0,z0);
    case 7 % a function u satisyinf uXn is not zero
        r = pi/4; x0 = -1; y0 = -1; z0 = -1;
        bm = 1; bp = 100; am = 1; ap = 200;
        pde = hyperIntf(am,ap,bm,bp,r,x0,y0,z0);
    case 8
        x0 = pi/5;
        bm = 1; bp = 10; am = 1; ap = 1;
        pde = ConstFun2(am,ap,bm,bp,x0);
    case 9 % twisted torus
        bm = 1; bp = 100; am = 1; ap = 200;
        domain = [-1.3,1.3,-1.3,1.3,-1.3,1.3];
        x1 = -0.3; y1 = 0; z1 = 0; r11 = pi/5; r12 = 0.2;
        x2 =  0.3; y2 = 0; z2 = 0; r21 = pi/5; r22 = 0.2;
        pde = elli3DtorusTwinHcurl2(am,ap,bm,bp,x1,y1,z1,r11,r12,x2,y2,z2,r21,r22);
      %%%%%%%%%%%%%%%%%%%% harmonic examples:
%     case 7 % harmonic
%         r1 = 1; r2 = 0; r3 = 0; x0 = pi/10^4; y0 = 0; z0 = 0;
%         bm = 1; bp = 10; am = 1; ap = 10;  kappa = 4;
%         pde = ConstFunHM(am,ap,bm,bp,kappa,r1,r2,r3,x0,y0,z0);
%         PDEoption.HM=1;
%     case 8  % ZJ's example （harmonic）
%         r1 = pi/5; r2 = pi/2;
%         n2 = 20; n1 = n2*(r2^2-r1^2);
%         bm = 1; bp = 10; am = 1; ap = 10; kappa = 4;
%         pde = elli3DcircIntf3HM(am,ap,bm,bp,kappa,r1,r2,n1,n2);
%         PDEoption.HM=1;
%     case 9
%         r = pi/5; x0 = 0; y0 = 0; z0 = 0;
%         bm = 1; bp = 10; am = 1; ap = 10; kappa = 4;
%         pde = LinearFun2HM(am,ap,bm,bp,kappa,r,x0,y0,z0);
%         PDEoption.HM=1;
end

%% Max Iteration
maxIt = 7;
time = zeros(9,maxIt);

%% three types of methods
PC = 'Bmg'; % use diagonal precondition
    
error = zeros(maxIt,3);
ratio = zeros(maxIt,3);
NumIte = zeros(maxIt,1);

for i = 1:maxIt
    
    %% 1. Generate Mesh and basic DOF
    tic
    nx = nx0 + 10*(i-1); h = (domain(2) - domain(1))/nx;
    ny = ny0 + 10*(i-1);
    nz = nz0 + 10*(i-1);
    %nx = NN(i); ny = nx; nz = nx; h = (domain(2) - domain(1))/nx;
    time(1,i) = nx;
    disp(' ')
    disp('*******************************************************************************')
    disp(['Partition =  ',int2str(nx),' X ',int2str(ny),' X ',int2str(nz)]);
    disp(' ')
    
    mesh = genMesh3D(domain, nx, ny, nz);
    mesh = enrichMesh3D(mesh,2); % Mesh detail level = 1 (for IFE).
    mesh = genIntfMesh3D(mesh,pde.intf);
    %%%%%%%%%% expansion   
    Explevel = 2;
    levelCount = 0;
    BasicEdgeNum = size(mesh.e,1);
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
    time(2,i) = toc;
    
    tic
    fem = genNedFEM3D(mesh,bc);
    disp(['number of interface element =  ', int2str(-min(mesh.tLoc)),...
        ', is ', num2str(100*-min(mesh.tLoc)/length(mesh.t)), '% of all elements']);
    time(3,i) = toc;
    
    %% 2. Generate FEMI data including IFE functions and their matrices
    tic
    meshI = genIVmesh(mesh);
    time(4,i) = toc;
    
    Voption.w1 = 1; %1 weight for the curl stabilization
    Voption.w2 = 1; %0.1 weight for the u stabilization
    femI =  genNedVIFEM3D(meshI,mesh,pde,Voption);
    time(5,i) = toc;
    
    
    %% 3. Assemble Matrix
    tic
    disp(' '); disp('Start Assembling Matrix');
    matrix = genMatCurlVIFE3D(pde,mesh,fem,meshI,femI,PDEoption);
    time(6,i) = toc;
    toc
    
    %% 4. Solve the linear system Au = f
    
    if strcmp(PC,'none')
        
        tic
        disp('Start Solving Linear System: using PCG with ichol precond');
        alpha = 10;
        L = ichol(matrix.A,struct('type','ict','droptol',1e-8,'diagcomp',alpha));
        [u,flag,relres,iter,resvec] = pcg(matrix.A,matrix.f,1e-9,2*10^4,L,L');
        disp(['number of iterations:', int2str(iter)]);
        time(7,i) = toc;
        
        tu = matrix.tu;
    elseif strcmp(PC,'BDPCG')
        
        tic
        disp('Start Solving Linear System: using PCG with block diagonal PCG');
        BasicEdgeNum = femI.BasicEdgeNum;
        D = diag(matrix.A(1:BasicEdgeNum,1:BasicEdgeNum));
        M = matrix.A(BasicEdgeNum+1:end,BasicEdgeNum+1:end);
        x0 = zeros(size(matrix.f));
        SolverMaxIter = 2*10^4;
        [L,U,P,Q] = lu(M);
        [x,info] = blockpcg(matrix.A,matrix.f,x0,1e-8,SolverMaxIter,D,L,U,P,Q);
        %[x,info] = blockpcg(matrix.A,matrix.f,x0,1e-8,SolverMaxIter,D,M);
        disp(['iteration number of PCG =',int2str(info.itStep)]);
        time(7,i) = toc;
        
        uh = x;
        tu = matrix.tu;
    elseif strcmp(PC,'Bmg')        
        %%%%%%%%%%%%%%
        % assemble the AP matrix: div(alpha grad) + I       
        tic
        femtype = 'P1';
        BC = [1,1,1,1,1,1]; % Dirichelet BC
        femH1 = genFEM3D(mesh,femtype,BC);
        opt.mass = 1; % generate the mass matrix
        opt.rhs = 0; % do not generate the risht hand-side
        opt.scl = 1;
        femIH1 =  genP1VIFEM3D(meshI,mesh,pde,opt);
        S = globMatrixVIFE3DStiff(pde.A,mesh,femH1,femH1);
        onefun = pde.B; %@(x,y,z) ones(size(x));
        M = globMatrixVIFE3DMass(onefun,mesh,femH1,femH1);
        if isfield(PDEoption,'HM')
            if PDEoption.HM ==1
                Atotal = S - kappa*M + femIH1.K + femIH1.S - kappa*femIH1.M;
            end
        else
            Atotal = S + M + femIH1.K + femIH1.S + femIH1.M;
        end
        node = meshI.node;
        NVdof = size(node,1);
        bcind = fem.bcind;
        [~,mapper] = boundaryNode3D(mesh,node,bcind);
        isBdNode = ones(NVdof,1);
        isBdNode(mapper) = 0;
        
        Tbd = spdiags(isBdNode,0,NVdof,NVdof);
        T = spdiags(1-isBdNode,0,NVdof,NVdof);
        AP = T*Atotal*T + Tbd;
        option.AP = AP;
        
        option.outsolver = 'cg';
        alpha = am*ones(size(femI.gdof,1),1);
        alpha(femI.eLoc>0) = ap;
        alpha(femI.eLoc==0) = (am+ap)/2;
        beta = bm*ones(size(femI.gdof,1),1);
        beta(femI.eLoc>0) = bp;
        beta(femI.eLoc==0) = (bm+bp)/2;
        option.alpha = alpha;
        option.beta = beta;
        option.solver = 'amg';
        edge = femI.gdof;
        option.isBdEdge = matrix.isBdEdge;
        option.smoother = 'BD';
        option.blklevel = 1;%ceil(sqrt(nx));
        BasicEdgeNum = femI.BasicEdgeNum;
        option.blkId = BasicEdgeNum;
        % size(femI.gdof,1): no submatrix for the direct solver at all
        % BasicEdgeNum: the basic submatrix for the direct solve
        option.fact = 'chol';
        %option.x0 = x0;
        
        [x,info] = amgMaxwellinterface(matrix.A,matrix.f,meshI.node,femI.gdof,option);
        uh = x;
        tu = matrix.tu;
        time(7,i) = toc;
        %uh = tu;
    end
    
    
    %% 4. Postprocess: Calculating Errors
    tic
    errND = max(abs(uh - tu)); % Error on nodes
    err.nd = errND; err.inf = 0; 
    err.l2 = 0; err.h1 = 0;
    if computErr == 1
        disp(' ');  disp('Start computing error in L2 norm');
        eNorm = 'L2'; disp(['Start computing error in ',eNorm,'  norm']);
        [errL2,errL2K1,errL2K2,errL2K3] = getCurlErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);
        
        eNorm = 'Curl'; disp(['Start computing error in ',eNorm,' norm']);
        [errCurl,errCurl1,errCurl2,errCurl3] = getCurlErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);
        err.l2 = errL2;
        err.curl = errCurl;
    end
    time(8,i) = toc;
    time(9,i) = 1e6*sum(time(2:6,i))/length(mesh.t);
    error(i,1) = errND; error(i,2) = errL2; error(i,3) = errCurl;
    NumIte(i) = info.itStep;
    
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
        ratio(i,1) = rNd;   ratio(i,2) = rL2; ratio(i,3) = rcurl;
        save(strcat('ratio_am',num2str(am),'ap',num2str(ap),'_bm',num2str(bm),'bp',num2str(bp)),'ratio')
        
    end
    err0 = err; h0 = h;
    save(strcat('error_am',num2str(am),'ap',num2str(ap),'_bm',num2str(bm),'bp',num2str(bp)),'error')
    
    disp(' '); disp('CPU Time')
    disp('   N     Mesh     FEM      MeshI   FEMI   Matrix   Solve    Error    Time/1M cell')
    formatSpec = '%4i  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f %7.2f  %7.2f\n';
    fprintf(formatSpec, time(:,1:i))
    
end

