function err = getErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm)
%% USAGE: calculate L^2 or semi-H^1 errors of FE/DG solution. 
%
% INPUTS:
% uh   --- a vector contains the FE solution or interpolation
% fun  --- exact function (u, Dx(u), Dy(u), or Dz(u))  
% mesh --- a struct data contains mesh information.
% fem  --- global degree of freedom for test function space
% type --- which norm 'L2','H1x','H1y','H1z'
%
% OUTPUTS:
% error --- error over the domain
%
% Last Modified: 09/18/2019 by Xu Zhang
% Last Modified: 07/02/2020 by Xu Zhang

%% 0. Initializaiton
if strcmp(fem.type,'P1') || strcmp(fem.type,'DGP1')||strcmp(fem.type,'CR') 
    feEvalBas = @evalP1Bas3D;
elseif strcmp(fem.type,'P2') || strcmp(fem.type,'DGP2')
    feEvalBas = @evalP2Bas3D;
end

if strcmp(eNorm,'L2')
    dind = [0,0,0];
    fun = pde.exactu;
    coefun = @(x,y,z) ones(size(x));
elseif strcmp(eNorm,'inf')
    dind = [0,0,0];
    fun = pde.exactu;
    coefun = @(x,y,z) ones(size(x));
elseif strcmp(eNorm,'H1x') || strcmp(eNorm,'Ex')
    dind = [1,0,0];
    fun = pde.Dxu;
    coefun = @(x,y,z) ones(size(x));
elseif strcmp(eNorm,'H1y') || strcmp(eNorm,'Ey')
    dind = [0,1,0];
    fun = pde.Dyu;
    coefun = @(x,y,z) ones(size(x));
elseif strcmp(eNorm,'H1z') || strcmp(eNorm,'Ez')
    dind = [0,0,1];
    fun = pde.Dzu;
    coefun = @(x,y,z) ones(size(x));
end

if strcmp(eNorm,'Ex') || strcmp(eNorm,'Ey') || strcmp(eNorm,'Ez')
    coefun = pde.A;
end

%% 1. Error on non-interface elements
dof = 4;  
ntID = find(mesh.tLoc > 0); 
AN = fem.area(ntID); gw = fem.gw;
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); 

uhK = uh(fem.t); uhKn = uhK(ntID,:); 
tuN = feval(fun, gxN, gyN, gzN);
coefN = feval(coefun, gxN, gyN, gzN);
uKN = 0; errK = zeros(size(mesh.t,1),1);
for i = 1:dof
    BAS = fem.bas(ntID,:,i);
    uKN = uKN + uhKn(:,i).*feEvalBas(BAS, gxN, gyN, gzN, dind);
end
errKN = sqrt(AN.*sum(coefN.^2.*(tuN-uKN).^2.*gw',2));
errK(ntID,1) = errKN;

%%% interface elements

NV = femI.NV;
i7 = find(NV==7);
i8 = find(NV==8);

node = meshI.node;
% transfer idx2cube to index of interface elements 1,2,3...
idx2cube = -mesh.tLoc(meshI.idx2cube);
tetElem = meshI.tetElem;
tetVolume = meshI.tetVolume;

if strcmp(eNorm,'inf')
    errI = max(max(abs(tuN-uKN)));
else
    errI = sum(errKN.^2);
end


for kk = 1:length(NV)
    
    nv = NV(kk);
    ElemID = femI.ElemID{kk}; TetID = ElemID(idx2cube);
    volume = tetVolume(TetID);
    ElemDoF = zeros(sum(TetID),nv);
    ElemDoF(ElemID,:) = femI.ElemDoF{kk}; 
    ElemDoF = ElemDoF(idx2cube(TetID),:);
    uhk = uh(ElemDoF);
    X1 = node(tetElem(TetID,1),:); 
    X2 = node(tetElem(TetID,2),:); 
    X3 = node(tetElem(TetID,3),:);
    X4 = node(tetElem(TetID,4),:);
    ng = 1;
    [gx, gy, gz] = gaussPtetra(X1, X2, X3, X4, ng);
    gw = gaussWtetra(ng);
    
    Bas1 = zeros(size(ElemID,1),3,nv); Bas2 = Bas1;
    Bas1(ElemID,:,:) = femI.Bas1{kk}; 
    Bas2(ElemID,:,:) = femI.Bas2{kk}; 
    Bas1 = Bas1(idx2cube(TetID),:,:); Bas2 = Bas2(idx2cube(TetID),:,:);
    Bas = zeros(sum(TetID),3,nv);
    %%%
    IFEconst = zeros(size(ElemID,1),nv);
    IFEconst(ElemID,:) = femI.IFEconst{kk};
    IFEconst = IFEconst(idx2cube(TetID),:);
    %%%
    Xm = zeros(size(ElemID,1),3);
    Xm(ElemID,:) = femI.Xmpt{kk};
    Xm = Xm(idx2cube(TetID),:);
    %%%
    piecetmp = meshI.tetElemLoc(TetID);
    Bas(piecetmp==1,:,:) = Bas1(piecetmp==1,:,:);
    Bas(piecetmp==2,:,:) = Bas2(piecetmp==2,:,:);
    uKI = zeros(size(gx));
    tuI = feval(fun, gx, gy, gz);
    
    if strcmp(eNorm,'L2')
        
        for i = 1:nv
            Basi = squeeze(Bas(:,:,i)); IFEconsti = IFEconst(:,i);
            uhi = (gx-Xm(:,1)).*Basi(:,1) + (gy-Xm(:,2)).*Basi(:,2) + (gz-Xm(:,3)).*Basi(:,3) + IFEconsti;
            uKI = uKI + uhk(:,i).*uhi;
        end
        errKI = sqrt(volume.*sum((tuI-uKI).^2.*gw',2));
        errI = errI + sum(errKI.^2);
        
    elseif strcmp(eNorm,'inf')
        
        for i = 1:nv
            Basi = squeeze(Bas(:,:,i)); IFEconsti = IFEconst(:,i);
            uhi = (gx-Xm(:,1)).*Basi(:,1) + (gy-Xm(:,2)).*Basi(:,2) + (gz-Xm(:,3)).*Basi(:,3) + IFEconsti;
            uKI = uKI + uhk(:,i).*uhi;
        end
        errKI = max(max(abs(tuI-uKI)));
        errI = max([errI,errKI]);
        
    else
        
        for i = 1:nv
            Basi = squeeze(Bas(:,:,i));
            uhi = Basi*dind';
            uKI = uKI + uhk(:,i).*uhi;
        end
        errKI = sqrt(volume.*sum((tuI-uKI).^2.*gw',2));
        errI = errI + sum(errKI.^2);
    end
    
end

if ~strcmp(eNorm,'inf')
    err = sqrt(errI);
else
    err = errI;
end

return