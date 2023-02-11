domain = [-1,1,-1,1,-1,1];
bc = [1,1,1,1,1,1]; % Dirichelet BC

% epsm = 8.854*10^(-12);
% epsp = 11.68*8.854*10^(-12);
% sigm = 10^5;
% sigp = 10^6;
% mum = 4*pi*10^(-7);
% mup = 4*pi*10^(-7);
% x0=0; y0=0; z0=-0.3; r1=pi/5; r2=0.2; CurrentT = 0; omega = 4;
% pde = TorusTime(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,CurrentT);

% epsm = 10;
% epsp = 1;
% sigm = 10;
% sigp = 1;
% mum = 10;
% mup = 1;
% x0=0; y0=0; z0=0; r1=pi/4; r2=pi/2; CurrentT = 0; omega = 2;
% n2 = 20; n1 = n2*(r2^2-r1^2);
% pde = elli3DcircIntf3Time(mum,mup,sigm,sigp,epsm,epsp,omega,r1,r2,n1,n2);
% TEND = 1;

epsm = 8.854*10^(-3);
epsp = 8.854*10^(-3);
sigm = 1;%100;
sigp = 1;%1;
mum = 4*pi;
mup = 4*pi;
x0=0; y0=0; z0=-0.3; r1=0.2; r2=pi/5; CurrentT = 0;
omega = 1; a = omega*sqrt(epsp*mup); b = 120; intPt = -1;
pde = TorusTimeInitial4(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,a,b,intPt);
TEND = 0.6;

nx0 = 10;
ny0 = nx0;
nz0 = nx0;
maxIt = 1;

for i = 1:maxIt
    
    CurrentT = 0;
    nx = nx0 + 10*(i-1); h = (domain(2) - domain(1))/nx;
    ny = ny0 + 10*(i-1);
    nz = nz0 + 10*(i-1);
    Ntime = ceil(TEND/sqrt(8.854*10^(-3)*4*pi))*10*nx;%ceil(sqrt(nx));
    deltaT = TEND/Ntime;
    
    mesh = genMesh3D(domain, nx, ny, nz);
    mesh = enrichMesh3D(mesh,2); % Mesh detail level = 1 (for IFE).
    mesh = genIntfMesh3D(mesh,pde.intf);
    
    Explevel = 3;
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
    
    fem = genNedFEM3D(mesh,bc);
    bm = epsm/deltaT^2 + sigm/(2*deltaT);
    bp = epsp/deltaT^2 + sigp/(2*deltaT);
    am = mum; ap = mup;
    femI = genNed1IFEM3D(mesh,fem,bm,bp,am,ap);
    
    S = globMatrixNedPGIFE3D(pde.Mu,1,mesh,femI,fem,fem);
    Me = globMatrixNedPGIFE3D(pde.Epslon,0,mesh,femI,fem,fem);
    Ms = globMatrixNedPGIFE3D(pde.Sig,0,mesh,femI,fem,fem);
    Atotal = Me/deltaT^2 + Ms/(2*deltaT) + S/2;
    
    Ndof = size(mesh.e,1);
    bdidx = zeros(Ndof,1);
    isBdEdge = true(Ndof,1);
    isBdEdge(fem.mapper) = false;
    bdidx(isBdEdge) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    A = T*Atotal*T + Tbd;
%     Me = T*Me*T + Tbd;
%     Ms = T*Ms*T + Tbd;
%     S = T*S*T + Tbd;
    
    exactu1Current = @(x,y,z) pde.exactu1(x,y,z,CurrentT);
    exactu2Current = @(x,y,z) pde.exactu2(x,y,z,CurrentT);
    exactu3Current = @(x,y,z) pde.exactu3(x,y,z,CurrentT);
    tu1 = sum(feval(exactu1Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
    tu2 = sum(feval(exactu2Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
    tu3 = sum(feval(exactu3Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
    % interface edges
    eID = find(mesh.eLoc<0);
    um1Current =  @(x,y,z) pde.um1(x,y,z,CurrentT);
    um2Current =  @(x,y,z) pde.um2(x,y,z,CurrentT);
    um3Current =  @(x,y,z) pde.um3(x,y,z,CurrentT);
    up1Current =  @(x,y,z) pde.up1(x,y,z,CurrentT);
    up2Current =  @(x,y,z) pde.up2(x,y,z,CurrentT);
    up3Current =  @(x,y,z) pde.up3(x,y,z,CurrentT);
    tu1I1 = sum(feval(um1Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
    tu2I1 = sum(feval(um2Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
    tu3I1 = sum(feval(um3Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
    tu1I2 = sum(feval(up1Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
    tu2I2 = sum(feval(up2Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
    tu3I2 = sum(feval(up3Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
    tu1(eID) = tu1I1 + tu1I2; tu2(eID) = tu2I1 + tu2I2; tu3(eID) = tu3I1 + tu3I2;
    tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
    tgt = tgt./sum(tgt.^2,2).^(1/2);
    uppre = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
    
    CurrentT = CurrentT + deltaT;
    exactu1Current = @(x,y,z) pde.exactu1(x,y,z,CurrentT);
    exactu2Current = @(x,y,z) pde.exactu2(x,y,z,CurrentT);
    exactu3Current = @(x,y,z) pde.exactu3(x,y,z,CurrentT);
    tu1 = sum(feval(exactu1Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
    tu2 = sum(feval(exactu2Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
    tu3 = sum(feval(exactu3Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
    % interface edges
    eID = find(mesh.eLoc<0);
    um1Current =  @(x,y,z) pde.um1(x,y,z,CurrentT);
    um2Current =  @(x,y,z) pde.um2(x,y,z,CurrentT);
    um3Current =  @(x,y,z) pde.um3(x,y,z,CurrentT);
    up1Current =  @(x,y,z) pde.up1(x,y,z,CurrentT);
    up2Current =  @(x,y,z) pde.up2(x,y,z,CurrentT);
    up3Current =  @(x,y,z) pde.up3(x,y,z,CurrentT);
    tu1I1 = sum(feval(um1Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
    tu2I1 = sum(feval(um2Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
    tu3I1 = sum(feval(um3Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
    tu1I2 = sum(feval(up1Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
    tu2I2 = sum(feval(up2Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
    tu3I2 = sum(feval(up3Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
    tu1(eID) = tu1I1 + tu1I2; tu2(eID) = tu2I1 + tu2I2; tu3(eID) = tu3I1 + tu3I2;
    tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
    tgt = tgt./sum(tgt.^2,2).^(1/2);
    upre = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
    
    while CurrentT < TEND
        
        CurrentT = CurrentT + deltaT;
        f1Current = @(x,y,z) pde.f1(x,y,z,CurrentT);
        f2Current = @(x,y,z) pde.f2(x,y,z,CurrentT);
        f3Current = @(x,y,z) pde.f3(x,y,z,CurrentT);
        rhsF1 = globNedRHSPGIFE3D(f1Current, mesh, fem, femI, 0, 1);
        rhsF2 = globNedRHSPGIFE3D(f2Current, mesh, fem, femI, 0, 2);
        rhsF3 = globNedRHSPGIFE3D(f3Current, mesh, fem, femI, 0, 3);
        rhsFCurrent = rhsF1 + rhsF2 + rhsF3;
        
        exactu1Current = @(x,y,z) pde.exactu1(x,y,z,CurrentT);
        exactu2Current = @(x,y,z) pde.exactu2(x,y,z,CurrentT);
        exactu3Current = @(x,y,z) pde.exactu3(x,y,z,CurrentT);
        tu1 = sum(feval(exactu1Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
        tu2 = sum(feval(exactu2Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
        tu3 = sum(feval(exactu3Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
        % interface edges
        eID = find(mesh.eLoc<0);
        um1Current =  @(x,y,z) pde.um1(x,y,z,CurrentT);
        um2Current =  @(x,y,z) pde.um2(x,y,z,CurrentT);
        um3Current =  @(x,y,z) pde.um3(x,y,z,CurrentT);
        up1Current =  @(x,y,z) pde.up1(x,y,z,CurrentT);
        up2Current =  @(x,y,z) pde.up2(x,y,z,CurrentT);
        up3Current =  @(x,y,z) pde.up3(x,y,z,CurrentT);
        tu1I1 = sum(feval(um1Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
        tu2I1 = sum(feval(um2Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
        tu3I1 = sum(feval(um3Current,femI.gex1,femI.gey1,femI.gez1).*femI.gew1,2);
        tu1I2 = sum(feval(up1Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
        tu2I2 = sum(feval(up2Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
        tu3I2 = sum(feval(up3Current,femI.gex2,femI.gey2,femI.gez2).*femI.gew2,2);
        tu1(eID) = tu1I1 + tu1I2; tu2(eID) = tu2I1 + tu2I2; tu3(eID) = tu3I1 + tu3I2;
        tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
        tgt = tgt./sum(tgt.^2,2).^(1/2);
        tuCurrent = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
        
        ub = tuCurrent;
        ub(fem.mapper) = 0;
        rhsB = Atotal*ub;
        JCurrent = rhsFCurrent - rhsB;
        %JCurrent(isBdEdge) = tuCurrent(isBdEdge);
        
        fcurrent = Me*(2*upre - uppre)/deltaT^2 + Ms*uppre/(2*deltaT) - S*uppre/2 + JCurrent;
        fcurrent(isBdEdge) = tuCurrent(isBdEdge);
        
        option.outsolver = 'gmres';
        Eplus = (mesh.eLoc==2);
        Eint = (mesh.eLoc<0);
        alpha = am*ones(size(mesh.e,1),1);
        alpha(Eplus) = ap;
        alpha(Eint) = (am+ap)/2;
        beta = bm*ones(size(mesh.e,1),1);
        beta(Eplus) = bp;
        beta(Eint) = (bm+bp)/2;
        option.alpha = alpha;
        option.beta = beta;
        option.solver = 'amg';
        edge = mesh.e;
        option.smoother = 'BD';
        option.blklevel = 0;
        option.blkId = BasicEdgeNum;
        [x,info] = amgMaxwellinterface2(A,fcurrent,mesh.p,edge,option);
        uppre = upre;
        upre = x;
        
    end
    
    uh = upre;
    
    disp(' ');  disp('Start computing error in L2 norm');
    eNorm = 'L2'; disp(['Start computing error in ',eNorm,'  norm']);
    [errL2,errL2K1,errL2K2,errL2K3] = getCurlErrIFE3DTime(uh, pde, mesh, fem, femI, eNorm, TEND);
    
    eNorm = 'Curl'; disp(['Start computing error in ',eNorm,' norm']);
    [errCurl,errCurl1,errCurl2,errCurl3] = getCurlErrIFE3DTime(uh, pde, mesh, fem, femI, eNorm, TEND);
    err.l2 = errL2;
    err.curl = errCurl;
    
    disp(' ')
    disp('Errors')
    disp('L2 norm     H1 norm')
    formatSpec = '%6.4e  %6.4e\n';
    fprintf(formatSpec, err.l2, err.curl)
    
    if i > 1
        format short
        rL2 = log(err0.l2/err.l2)./log(h0/h);
        rcurl = log(err0.curl/err.curl)./log(h0/h);
        disp(' ')
        disp('Convergence Rate')
        disp('L2 norm     Curl norm')
        formatSpec = '%6.4f      %6.4f\n';
        fprintf(formatSpec, rL2, rcurl)
    end
    err0 = err; h0 = h;
    
end