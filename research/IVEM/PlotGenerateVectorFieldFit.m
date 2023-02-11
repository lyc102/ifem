%load('BoxTorusMesh.mat')
%meshFit = mesh;
uh = UHfit(:,22);
meshFit = enrichMesh3D(meshFit,0);
bc = [1,1,1,1,1,1];
femFit = genNedFEM3D(meshFit,bc);
ng = 1;
X1 = meshFit.p(meshFit.t(:,1),:); X2 = meshFit.p(meshFit.t(:,2),:); 
X3 = meshFit.p(meshFit.t(:,3),:); X4 = meshFit.p(meshFit.t(:,4),:); 
gw = gaussWtetra(ng);
[gx,gy,gz] = gaussPtetra(X1,X2,X3,X4,ng);
feEvalBas = @EvalNed1Bas3D;

uhK = uh(femFit.g2ldof); 
uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:6
    uK1 = uK1 + uhK(:,i).*feEvalBas(femFit.bas, ':', gx, gy, gz, i, 0, 1).*femFit.t_e_orit(:,i);
    uK2 = uK2 + uhK(:,i).*feEvalBas(femFit.bas, ':', gx, gy, gz, i, 0, 2).*femFit.t_e_orit(:,i);
    uK3 = uK3 + uhK(:,i).*feEvalBas(femFit.bas, ':', gx, gy, gz, i, 0, 3).*femFit.t_e_orit(:,i);
end

Vhfit = [uK1,uK2,uK3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% separate evaluation
domain = [-1,1,-1,1,-1,1];
nx = 64;
ny = nx;
nz = nx;
mesh = genMesh3D(domain, nx, ny, nz);
mesh = enrichMesh3D(mesh,2);
tol = 10^(-12);
zCut = -0.3125;
fp1 = mesh.p(mesh.f(:,1),:);
fp2 = mesh.p(mesh.f(:,2),:);
fp3 = mesh.p(mesh.f(:,3),:);
fZid = (abs(fp1(:,3)-zCut)<tol).*(abs(fp2(:,3)-zCut)<tol).*(abs(fp3(:,3)-zCut)<tol);
fZid = find(fZid==1);
fmpt = (fp1(fZid,:)+fp2(fZid,:)+fp3(fZid,:))/3;

fMallfit = (meshFit.p(meshFit.t(:,1),:)+meshFit.p(meshFit.t(:,2),:)+...
    meshFit.p(meshFit.t(:,3),:)+meshFit.p(meshFit.t(:,4),:))/4;

tic
[tr,tj]=findtria(meshFit.p,meshFit.t,fmpt);
toc

VhCut = [uK1(tj),uK2(tj),uK3(tj)];

fmptPlan = fmpt(:,1:2);
zDT = delaunay(fmptPlan);


%%%%% trisulf plot
trisurf(zDT,fmptPlan(:,1),fmptPlan(:,2),VhCut(:,1));
shading interp

%%%%% vector field plot
quiver3(fmpt(:,1),fmpt(:,2),fmpt(:,3),VhCut(:,1),VhCut(:,2),VhCut(:,3))

%%%%%% contour plot
tricontour(fmptPlan,zDT,VhCut(:,1),10)

tid1 = find(mesh.tLoc==1);
tid1 = tid1(1:50:end);
fmpt1 = fMall(tid1,:);
[tr,tj]=findtria(meshFit.p,meshFit.t,fmpt1);
Vh1 = Vhfit(tj,:);
quiver5(fmpt1(:,1),fmpt1(:,2),fmpt1(:,3),...
    Vh1(:,1),Vh1(:,2),Vh1(:,3),5,'filled');
