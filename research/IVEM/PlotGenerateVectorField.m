Vh = zeros(size(mesh.t,1),3);

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
uhK = uh(fem.g2ldof); 
uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:dof
    uK1 = uK1 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, 0, 1).*fem.t_e_orit(ntID,i);
    uK2 = uK2 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, 0, 2).*fem.t_e_orit(ntID,i);
    uK3 = uK3 + uhK(ntID,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, 0, 3).*fem.t_e_orit(ntID,i);
end

Vh(ntID,:) = [uK1,uK2,uK3];

% interface elements
epsm = 8.85*10^(-3)*5;
epsp = 8.85*10^(-3);
sigm = 100;
sigp = 1;
mum = (4*pi);
mup = (4*pi)*2;
bm = epsm/deltaT^2 + sigm/(2*deltaT);
bp = epsp/deltaT^2 + sigp/(2*deltaT);
am = mum^(-1); ap = mup^(-1);
femI = genNed1IFEM3D(mesh,fem,am,ap,bm,bp,option);
gx = femI.gx; gy = femI.gy; gz = femI.gz; 
tIntfID = femI.tIntfID;
itID = find(mesh.tLoc < 0); ntI = length(itID);
uhK = uh(fem.g2ldof); uhKitmp = uhK(itID,:); uhKi = zeros(length(femI.g2ldof),6);
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
Vh(itID,:) = [Vi1,Vi2,Vi3];

fMall = (mesh.p(mesh.t(:,1),:)+mesh.p(mesh.t(:,2),:)+...
    mesh.p(mesh.t(:,3),:)+mesh.p(mesh.t(:,4),:))/4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 10^(-12);
zCut = -0.3125;
fp1 = mesh.p(mesh.f(:,1),:);
fp2 = mesh.p(mesh.f(:,2),:);
fp3 = mesh.p(mesh.f(:,3),:);
fZid = (abs(fp1(:,3)-zCut)<tol).*(abs(fp2(:,3)-zCut)<tol).*(abs(fp3(:,3)-zCut)<tol);
fZid = find(fZid==1);
fmpt = (fp1(fZid,:)+fp2(fZid,:)+fp3(fZid,:))/3;
fmptPlan = fmpt(:,1:2);
zDT = delaunay(fmptPlan);

fZTid = mesh.f_t(fZid,1);
VhCut = Vh(fZTid,:);

%%%%% trisulf plot
trisurf(zDT,fmptPlan(:,1),fmptPlan(:,2),VhCut(:,1));
shading interp

%%%%%% contour plot
tricontour(fmptPlan,zDT,VhCut(:,1),10)

%%%%% vector field plot
quiver3(fmpt(:,1),fmpt(:,2),fmpt(:,3),VhCut(:,1),VhCut(:,2),VhCut(:,3))

tid1 = find(mesh.tLoc==1);
tid1 = tid1(1:50:end);
quiver3(fMall(tid1,1),fMall(tid1,2),fMall(tid1,3),Vh(tid1,1),Vh(tid1,2),Vh(tid1,3),20)

tid1 = find(mesh.tLoc==1);
tid1 = tid1(1:50:end);
quiver5(fMall(tid1,1),fMall(tid1,2),fMall(tid1,3),...
    Vh(tid1,1),Vh(tid1,2),Vh(tid1,3),5);

tidzn = find(fMall(:,3)<0);
tidzn = tidzn(1:100:end);
quiver5(fMall(tidzn,1),fMall(tidzn,2),fMall(tidzn,3),...
    Vh(tidzn,1),Vh(tidzn,2),Vh(tidzn,3),3);

