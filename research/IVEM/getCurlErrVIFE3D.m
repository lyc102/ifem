function [err,errK1,errK2,errK3] = getCurlErrVIFE3D(uh, pde, mesh, meshI, femI, fem, eNorm);

%% USAGE: calculate L^2 or semi-H^1 errors of FE/DG solution. 
%
% INPUTS:
% uh   --- a vector contains the FE solution or interpolation
% fun  --- exact function (u, Dx(u), Dy(u), or Dz(u))  
% fem  --- global degree of freedom for test function space
% type --- which norm 'L2','H1x','H1y','H1z'
%
% OUTPUTS:
% error --- error over the domain
%
% Last Modified: 08/07/2020 by Xu Zhang

feEvalBas = @EvalNed1Bas3D;

%% 1. Error on non-interface elements 
dof = 6;  
ntID = find(mesh.tLoc>0); 
AN = fem.area(ntID); gw = fem.gw;
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:);
errK1 = zeros(size(mesh.t,1),1); errK2 = errK1; errK3 = errK1;

if strcmp(eNorm,'L2')
    dind = 0;
    fun1 = pde.exactu1; fun2 = pde.exactu2; fun3 = pde.exactu3;
elseif strcmp(eNorm,'Curl')
    dind = 1;
    fun1 = pde.Dxu; fun2 = pde.Dyu; fun3 = pde.Dzu;
end
tu1 = feval(fun1, gxN, gyN, gzN);
tu2 = feval(fun2, gxN, gyN, gzN);
tu3 = feval(fun3, gxN, gyN, gzN);
uhK = uh(femI.g2ldofNint); 

uK1 = 0; uK2 = 0; uK3 = 0;
for i = 1:dof
    uK1 = uK1 + uhK(:,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, 1).*fem.t_e_orit(ntID,i);
    uK2 = uK2 + uhK(:,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, 2).*fem.t_e_orit(ntID,i);
    uK3 = uK3 + uhK(:,i).*feEvalBas(fem.bas, ntID, gxN, gyN, gzN, i, dind, 3).*fem.t_e_orit(ntID,i);
end
errK1(ntID) = sqrt(AN.*sum((tu1-uK1).^2.*gw',2));
errK2(ntID) = sqrt(AN.*sum((tu2-uK2).^2.*gw',2));
errK3(ntID) = sqrt(AN.*sum((tu3-uK3).^2.*gw',2));

errN = sum((errK1.^2+errK2.^2+errK3.^2)).^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Error on interface elements 
NE = femI.une;
node = meshI.node;
% transfer idx2cube to index of interface elements 1,2,3...
idx2cube = -mesh.tLoc(meshI.idx2cube);
tetElem = meshI.tetElem;
tetVolume = meshI.tetVolume;
ntI = sum(mesh.tLoc<0); 

errI = sum(errN.^2);


for kk = 1:length(NE)
    
    ne = NE(kk);
    ElemID = femI.ElemID{kk}; TetID = ElemID(idx2cube);
    volume = tetVolume(TetID);
    ElemDoF = zeros(sum(TetID),ne);
    ElemDoF(ElemID,:) = femI.ElemDoF{kk}; 
    ElemDoF = ElemDoF(idx2cube(TetID),:);
    uhk = uh(ElemDoF);
    X1 = node(tetElem(TetID,1),:); 
    X2 = node(tetElem(TetID,2),:); 
    X3 = node(tetElem(TetID,3),:);
    X4 = node(tetElem(TetID,4),:);
    ng = 4;
    [gx, gy, gz] = gaussPtetra(X1, X2, X3, X4, ng);
    gw = gaussWtetra(ng);
    
    Bas1 = zeros(size(ElemID,1),3,ne); Bas2 = Bas1;
    if strcmp(eNorm,'L2')
        Bas1(ElemID,:,:) = femI.BasB1{kk};
        Bas2(ElemID,:,:) = femI.BasB2{kk};
    elseif strcmp(eNorm,'Curl')
        Bas1(ElemID,:,:) = femI.BasA1{kk};
        Bas2(ElemID,:,:) = femI.BasA2{kk};
    end
    Bas1 = Bas1(idx2cube(TetID),:,:); Bas2 = Bas2(idx2cube(TetID),:,:);
    Bas = zeros(sum(TetID),3,ne);
    %%%
    piecetmp = meshI.tetElemLoc(TetID);
    Bas(piecetmp==1,:,:) = Bas1(piecetmp==1,:,:);
    Bas(piecetmp==2,:,:) = Bas2(piecetmp==2,:,:);
    uKI = zeros(size(gx,1),3);
    tu1I1 = feval(fun1, gx, gy, gz);
    tu1I2 = feval(fun2, gx, gy, gz);
    tu1I3 = feval(fun3, gx, gy, gz);
    
    
    for i = 1:ne
        Basi = squeeze(Bas(:,:,i));
        uKI = uKI + uhk(:,i).*Basi;
    end
    errKI = sqrt(volume.*sum( (tu1I1-uKI(:,1)).^2.*gw'+...
        (tu1I2-uKI(:,2)).^2.*gw' + (tu1I3-uKI(:,3)).^2.*gw',2));
    errI = errI + sum(errKI.^2);
    
end

if ~strcmp(eNorm,'inf')
    err = sqrt(errI);%*ntI^(1/4);
else
    err = errI;
end

return