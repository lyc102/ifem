domain = [-1,1,-1,1,-1,1];
bc = [1,1,1,1,1,1]; % Dirichelet BC

load('BoxTorusMesh.mat')
% epsm = 8.85*10^(-5)*5;
% epsp = 8.85*10^(-5);
% sigm = 100;
% sigp = 1;
% mum = (4*pi);
% mup = (4*pi)*2;
% x0=0; y0=0; z0=-0.3; r1=0.2; r2=pi/5; omega = 20; stre = 1;
% pde = TorusTime1(mum,mup,sigm,sigp,epsm,epsp,omega,stre,x0,y0,z0,r1,r2);
% TEND = 0.1;
epsm = 8.854*10^(-2);
epsp = 8.854*10^(-3);
sigm = 10;%100;
sigp = 1;%1;
mum = 4*pi*5;
mup = 4*pi;
x0=0; y0=0; z0=-0.3; r1=0.2; r2=pi/5; omega = 1; a = omega*sqrt(epsp*mup); b= 120; intPt = -1;
pde = TorusTimeInitial3(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,a,b,intPt);
TEND = 0.8;
InitCond = 1;
BDCond = 1;

nx = 64;  h=(domain(2) - domain(1))/nx;
ny = nx;
nz = nx;
    
Ntime = ceil(TEND/sqrt(8.854*10^(-3)*4*pi*2))*5*nx;%ceil(sqrt(nx));
deltaT = TEND/Ntime;
mesh = enrichMesh3D(mesh,0);
fem = genNedFEM3D(mesh,bc);
bm = epsm/deltaT^2 + sigm/(2*deltaT);
bp = epsp/deltaT^2 + sigp/(2*deltaT);
am = mum^(-1); ap = mup^(-1);

S = globMatrixNedFit3D(am,ap,1,mesh,fem,fem);
Me = globMatrixNedFit3D(epsm,epsp,0,mesh,fem,fem);
Ms = globMatrixNedFit3D(sigm,sigp,0,mesh,fem,fem);
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

UHfit =zeros(size(mesh.e,1),ceil(Ntime/10)+2);
NumStep = 0;
CurrentT = 0;
SaveCount = 0;

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
NEdof = size(A,1);

while CurrentT < TEND
    
    NumStep = NumStep + 1;
    CurrentT = CurrentT + deltaT;
    fm1Current = @(x,y,z) pde.fm1(x,y,z,CurrentT);
    fm2Current = @(x,y,z) pde.fm2(x,y,z,CurrentT);
    fm3Current = @(x,y,z) pde.fm3(x,y,z,CurrentT);
    fp1Current = @(x,y,z) pde.fp1(x,y,z,CurrentT);
    fp2Current = @(x,y,z) pde.fp2(x,y,z,CurrentT);
    fp3Current = @(x,y,z) pde.fp3(x,y,z,CurrentT);
    rhsF1 = globNedFitRHS3D(fm1Current,fp1Current, mesh, fem, 0, 1);
    rhsF2 = globNedFitRHS3D(fm2Current,fp2Current, mesh, fem, 0, 2);
    rhsF3 = globNedFitRHS3D(fm3Current,fp3Current, mesh, fem, 0, 3);
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
    end
    
    ub = tuCurrent;
    ub(fem.mapper) = 0;
    rhsB = Atotal*ub;
    JCurrent = rhsFCurrent - rhsB;
    %JCurrent(isBdEdge) = tuCurrent(isBdEdge);
    
    fcurrent = Me*(2*upre - uppre)/deltaT^2 + Ms*uppre/(2*deltaT) - S*uppre/2 + JCurrent;
    fcurrent(isBdEdge) = tuCurrent(isBdEdge);
    
    option.outsolver = 'cg';
    option.alpha = alpha;
    option.beta = beta;
    option.solver = 'amg';
    edge = mesh.e;
    option.smoother = 'BD';
    option.blklevel = 0;
    option.blkId = NEdof;
    [x,info] = amgMaxwellinterface2(A,fcurrent,mesh.p,edge,option);
    uppre = upre;
    upre = x;
    
    if mod(NumStep,40) == 0
        UHfit(:,SaveCount+1) = uppre;
        UHfit(:,SaveCount+2) = upre;
        SaveCount = SaveCount+2;
        save('UHfit1','UHfit','-v7.3');
        disp('save data')
    end
    
end

save('UHfit1','UHfit'); 