domain = [-1,1,-1,1,-1,1];
bc = [1,1,1,1,1,1]; % Dirichelet BC

% epsm = 8.854*10^(-12);
% epsp = 11.68*8.854*10^(-12);
% sigm = 10^5;
% sigp = 10^6;
% mum = 4*pi*10^(-7);
% mup = 4*pi*10^(-7);
epsm = 8.85*10^(-2)*5;
epsp = 8.85*10^(-2);
sigm = 100;
sigp = 1;
mum = (4*pi)^(-1);
mup = (4*pi)^(-1)*2;
x0=0; y0=0; z0=-0.3; r1=0.2; r2=pi/5; omega = 20; a = 1; b= 150; intPt = -1;
pde = TorusTimeInitial2(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,a,b,intPt);
TEND = 0.1;
InitCond = 1;
BDCond = 1;

nx = 64;  h=(domain(2) - domain(1))/nx;
ny = nx;
nz = nx;
    
Ntime = 5*nx;%ceil(sqrt(nx));
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
am = mum^(-1); ap = mup^(-1);
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

if InitCond == 0
    uppre = zeros(size(mesh.e,1),1);
    upre = zeros(size(mesh.e,1),1);
elseif InitCond ~=0
    tu1 = sum(feval(pde.E1,fem.gex,fem.gey,fem.gez,0).*fem.gew,2);
    tu2 = sum(feval(pde.E2,fem.gex,fem.gey,fem.gez,0).*fem.gew,2);
    tu3 = sum(feval(pde.E3,fem.gex,fem.gey,fem.gez,0).*fem.gew,2);
    tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
    tgt = tgt./sum(tgt.^2,2).^(1/2);
    uppre = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
    tu1 = sum(feval(pde.Et1,fem.gex,fem.gey,fem.gez,0).*fem.gew,2);
    tu2 = sum(feval(pde.Et2,fem.gex,fem.gey,fem.gez,0).*fem.gew,2);
    tu3 = sum(feval(pde.Et3,fem.gex,fem.gey,fem.gez,0).*fem.gew,2);
    upretmp = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
    upre = uppre + deltaT*upretmp;
end
UH =zeros(size(mesh.e,1),ceil(Ntime/10)+2);
NumStep = 0;
CurrentT = 0;
SaveCount = 0;
pmt = (mesh.p(mesh.e(:,2),:)+mesh.p(mesh.e(:,1),:))/2;
pf0 = (abs(pmt(:,1)-domain(2))<10^(-10));

while CurrentT < TEND
    
    NumStep = NumStep + 1;
    CurrentT = CurrentT + deltaT;
    f1Current = @(x,y,z) pde.f1(x,y,z,CurrentT);
    f2Current = @(x,y,z) pde.f2(x,y,z,CurrentT);
    f3Current = @(x,y,z) pde.f3(x,y,z,CurrentT);
    rhsF1 = globNedRHSPGIFE3D(f1Current, mesh, fem, femI, 0, 1);
    rhsF2 = globNedRHSPGIFE3D(f2Current, mesh, fem, femI, 0, 2);
    rhsF3 = globNedRHSPGIFE3D(f3Current, mesh, fem, femI, 0, 3);
    rhsFCurrent = rhsF1 + rhsF2 + rhsF3;
    
    if BDCond == 0
        tuCurrent = zeros(size(mesh.e,1),1);
    elseif BDCond ~= 0
        exactu1Current = @(x,y,z) pde.E1(x,y,z,CurrentT);
        exactu2Current = @(x,y,z) pde.E2(x,y,z,CurrentT);
        exactu3Current = @(x,y,z) pde.E3(x,y,z,CurrentT);
        tu1 = sum(feval(exactu1Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
        tu2 = sum(feval(exactu2Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
        tu3 = sum(feval(exactu3Current,fem.gex,fem.gey,fem.gez).*fem.gew,2);
        tgt = mesh.p(mesh.e(:,2),:) - mesh.p(mesh.e(:,1),:);
        tgt = tgt./sum(tgt.^2,2).^(1/2);
        tuCurrent = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
        pmt(pf0) = 0;
    end
    
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
    
    if mod(NumStep,10) == 0
       SaveCount = SaveCount+1;
       UH(:,SaveCount) = x;
       disp('save data')
       save('UH5','UH'); 
    end
    
end

save('UH5','UH'); 
    
%     uh = upre;
%     
%     disp(' ');  disp('Start computing error in L2 norm');
%     eNorm = 'L2'; disp(['Start computing error in ',eNorm,'  norm']);
%     [errL2,errL2K1,errL2K2,errL2K3] = getCurlErrIFE3DTime(uh, pde, mesh, fem, femI, eNorm, TEND);
%     
%     eNorm = 'Curl'; disp(['Start computing error in ',eNorm,' norm']);
%     [errCurl,errCurl1,errCurl2,errCurl3] = getCurlErrIFE3DTime(uh, pde, mesh, fem, femI, eNorm, TEND);
%     err.l2 = errL2;
%     err.curl = errCurl;
%     
%     disp(' ')
%     disp('Errors')
%     disp('L2 norm     H1 norm')
%     formatSpec = '%6.4e  %6.4e\n';
%     fprintf(formatSpec, err.l2, err.curl)
%     
%     if i > 1
%         format short
%         rL2 = log(err0.l2/err.l2)./log(h0/h);
%         rcurl = log(err0.curl/err.curl)./log(h0/h);
%         disp(' ')
%         disp('Convergence Rate')
%         disp('L2 norm     Curl norm')
%         formatSpec = '%6.4f      %6.4f\n';
%         fprintf(formatSpec, rL2, rcurl)
%     end
%     err0 = err; h0 = h;
    