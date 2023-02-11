if exist('meshFit','var') == 0
    load('BoxTorusMesh.mat')
    meshFit = mesh;
    meshFit = enrichMesh3D(meshFit,0);
end
step = 24;
uhIFE  = UH(:,step);
if exist('UHfit','var') == 0
    uhFE = zeros(size(meshFit.e,1),1);
else
    uhFE = UHfit(:,step);
end

%%%%%%% setup
domain = [-1,1,-1,1,-1,1];
bc = [1,1,1,1,1,1]; 
nx = 64;  h=(domain(2) - domain(1))/nx;
ny = nx;
nz = nx;

epsm = 5*10^(-2)*2;
epsp = 5*10^(-2);
sigm = 1;
sigp = 0.1;
mum = (4*pi)*3;
mup = (4*pi);
x0=0; y0=0; z0=-0.3; r1=0.2; r2=pi/5; omega = 1; a = omega*sqrt(epsp*mup); b= 150; intPt = -1;
pde = TorusTimeInitial3(mum,mup,sigm,sigp,epsm,epsp,omega,x0,y0,z0,r1,r2,a,b,intPt);
TEND = 1.5;
Ntime = ceil(TEND/sqrt(mup*epsp))*nx;%ceil(sqrt(nx));
deltaT = TEND/Ntime;

%%%%%%%%%%%%%%%%%%%%% for IFE functions

mesh = genMesh3D(domain, nx, ny, nz);
mesh = enrichMesh3D(mesh,2); % Mesh detail level = 1 (for IFE).
mesh = genIntfMesh3D(mesh,pde.intf);

% level 1
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

feEvalBas = @EvalNed1Bas3D;
IfeEvalBas = @evalNed1IFEBas3D;
option.ng = 1;
option.nge = 1;

p = mesh.p; t = mesh.t;
fem = genNedFEM3D(mesh,bc,option);
ntID = find(mesh.tLoc>0); 
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:);

dof = fem.ldof;  
ntID = find(mesh.tLoc>0); 
uhK = uhIFE(fem.g2ldof); 
uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:dof
    uK1 = uK1 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, 0, 1).*fem.t_e_orit(ntID,i);
    uK2 = uK2 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, 0, 2).*fem.t_e_orit(ntID,i);
    uK3 = uK3 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, 0, 3).*fem.t_e_orit(ntID,i);
end

nintID = length(ntID);
VhIFE(ntID,:) = [uK1,uK2,uK3];
%VhIFE(1:nintID,:) = [uK1,uK2,uK3];

bm = epsm/deltaT^2 + sigm/(2*deltaT);
bp = epsp/deltaT^2 + sigp/(2*deltaT);
am = mum^(-1); ap = mup^(-1);
%am = 1; ap = 1;
femI = genNed1IFEM3D(mesh,fem,bm,bp,am,ap,option);
gx = femI.gx; gy = femI.gy; gz = femI.gz; 
tIntfID = femI.tIntfID;
itID = find(mesh.tLoc < 0); ntI = length(itID);
uhK = uhIFE(fem.g2ldof); uhKitmp = uhK(itID,:); uhKi = zeros(length(femI.g2ldof),6);
uK1 = 0; uK2 = 0; uK3 = 0;
id = 0;
for i = 1:ntI
    tmp = uhKitmp(i,:);
    uhKi(id+1:id+femI.quadT(i),:) = repmat(tmp,femI.quadT(i),1);
    id = id + femI.quadT(i);
end

for i = 1:dof
    BAS = femI.bas(:,:,i);
    uK1 = uK1 + uhKi(:,i).*IfeEvalBas(BAS, gx, gy, gz, 0, 1).*fem.t_e_orit(tIntfID,i);
    uK2 = uK2 + uhKi(:,i).*IfeEvalBas(BAS, gx, gy, gz, 0, 2).*fem.t_e_orit(tIntfID,i);
    uK3 = uK3 + uhKi(:,i).*IfeEvalBas(BAS, gx, gy, gz, 0, 3).*fem.t_e_orit(tIntfID,i);
end
tIntfIDtmp = -mesh.tLoc(tIntfID);
Vi1 = accumarray(tIntfIDtmp,uK1)./femI.quadT;
Vi2 = accumarray(tIntfIDtmp,uK2)./femI.quadT;
Vi3 = accumarray(tIntfIDtmp,uK3)./femI.quadT;
% qq = femI.quadT; qq(1)=1; qq =cumsum(qq);
% Vi1 = uK1(qq); Vi2 = uK2(qq); Vi3 = uK3(qq);
VhIFE(itID,:) = [Vi1,Vi2,Vi3];
%VhIFE(nintID+1:nintID+size(uK1,1),:) = [uK1,uK2,uK3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for FE functions

bc = [1,1,1,1,1,1];
femFit = genNedFEM3D(meshFit,bc);
ng = 1;
X1 = meshFit.p(meshFit.t(:,1),:); X2 = meshFit.p(meshFit.t(:,2),:); 
X3 = meshFit.p(meshFit.t(:,3),:); X4 = meshFit.p(meshFit.t(:,4),:); 
gw = gaussWtetra(ng);
[gx,gy,gz] = gaussPtetra(X1,X2,X3,X4,ng);
feEvalBas = @EvalNed1Bas3D;

uhK = uhFE(femFit.g2ldof); 
uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:6
    uK1 = uK1 + uhK(:,i).*feEvalBas(femFit.bas, ':', gx, gy, gz, i, 0, 1).*femFit.t_e_orit(:,i);
    uK2 = uK2 + uhK(:,i).*feEvalBas(femFit.bas, ':', gx, gy, gz, i, 0, 2).*femFit.t_e_orit(:,i);
    uK3 = uK3 + uhK(:,i).*feEvalBas(femFit.bas, ':', gx, gy, gz, i, 0, 3).*femFit.t_e_orit(:,i);
end

VhFE = [uK1,uK2,uK3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot

corInd = 3;
tol = 10^(-12);
zCut = -0.3125;
fp1 = mesh.p(mesh.f(:,1),:);
fp2 = mesh.p(mesh.f(:,2),:);
fp3 = mesh.p(mesh.f(:,3),:);
fZid = (abs(fp1(:,corInd)-zCut)<tol).*(abs(fp2(:,corInd)-zCut)<tol).*(abs(fp3(:,corInd)-zCut)<tol);
fZid = find(fZid==1);
fmpt = (fp1(fZid,:)+fp2(fZid,:)+fp3(fZid,:))/3;
fmptPlan = fmpt(:,setdiff([1,2,3],corInd));
zDT = delaunay(fmptPlan);

fZTid = mesh.f_t(fZid,1);
VhIFECut = VhIFE(fZTid,:);

[tr,tj]=findtria(meshFit.p,meshFit.t,fmpt);
VhFECut = VhFE(tj,:);

%%%%%%%%%%%%% contour plot
dsty = 16;
figure(1)
tricontour(fmptPlan,zDT,VhIFECut(:,1),dsty)
figure(2)
tricontour(fmptPlan,zDT,VhFECut(:,1),dsty)

%%%%%%%%%%%% trisurf plot
figure(1)
trisurf(zDT,fmptPlan(:,1),fmptPlan(:,2),VhIFECut(:,2));
shading interp
colormap('jet')
caxis([min(VhFECut(:,2)) max(VhFECut(:,2))])
figure(2)
trisurf(zDT,fmptPlan(:,1),fmptPlan(:,2),VhFECut(:,2));
shading interp
colormap('jet')
caxis([min(VhFECut(:,2)) max(VhFECut(:,2))])


%%%%%%%%%%%% vector field plot
fMall = (mesh.p(mesh.t(:,1),:)+mesh.p(mesh.t(:,2),:)+...
    mesh.p(mesh.t(:,3),:)+mesh.p(mesh.t(:,4),:))/4;
tid1 = find(mesh.tLoc==1);
tid1 = tid1(1:50:end);
fmpt1 = fMall(tid1,:);
Vh1IFE = VhIFE(tid1,:);

[tr,tj]=findtria(meshFit.p,meshFit.t,fmpt1);
Vh1FE = VhFE(tj,:);

mgnt1 = (Vh1IFE(:,1).^2+Vh1IFE(:,2).^2+Vh1IFE(:,3).^2).^(1/2);
mgnt2 = (Vh1FE(:,1).^2+Vh1FE(:,2).^2+Vh1FE(:,3).^2).^(1/2);
mgnmax = max([max(mgnt1),max(mgnt2)]);
mgnt1 = mgnt1/mgnmax;
mgnt2 = mgnt2/mgnmax;

figure(1)
sz = 5;
quiver3d(fmpt1(:,1),fmpt1(:,2),fmpt1(:,3),...
    Vh1IFE(:,1),Vh1IFE(:,2),Vh1IFE(:,3),mgnt1,sz);
axis([-1, 1, -1, 1, -1, 0])
view(-20,30)
box on
%colormap('jet');

figure(2)
quiver3d(fmpt1(:,1),fmpt1(:,2),fmpt1(:,3),...
    Vh1FE(:,1),Vh1FE(:,2),Vh1FE(:,3),mgnt2,sz);
axis([-1, 1, -1, 1, -1, 0])
view(-20,30)
box on
%colormap('jet');




