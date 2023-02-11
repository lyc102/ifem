function femI = genP1VIFEM3D(meshI,mesh,pde,option)
%% use face and face2elem data structure to compute ife functions
bm = pde.bm; bp = pde.bp;

time = 0;
tic 
vSign = [mesh.pLoc;zeros(size(mesh.eIntP,1),1)];
node = [mesh.p;mesh.eIntP];
face = meshI.tface;
face2elem = meshI.tface2elem;
N = size(node, 1); Ndof = N;
NC = size(mesh.t, 1);

NP = max(face2elem);
isCutPoly = false(NP, 1);
isCutPoly(face2elem) = true;
% this is the same as:
% isCutPoly = (mesh.tLoc<0);
% isCutPoly = isCutPoly(1:max(face2elem));
cutPolyIdx = zeros(NP, 1);
isExtPoly = false(NP,1);
isExtPoly(NC+1:NP) = true;

NP = sum(isCutPoly); % this new NP is max(-mesh.tLoc(mesh.tLoc<0));
cutPolyIdx(isCutPoly) = 1:NP;
% this is the same as:
% cutPolyIdx = zeros(size(mesh.tLoc,1), 1);
% cutPolyIdx(mesh.tLoc<0) = -mesh.tLoc(mesh.tLoc<0);
% cutPolyIdx = cutPolyIdxNew(1:max(face2elem));
face2elem = cutPolyIdx(face2elem);
% the old face2elem is the global index of interface elements
% the new one is only the index of interface elements
isExtPoly = isExtPoly(isCutPoly);
%%%clear isCutPoly cutPolyIdx

% prepare new data structure to compute projections
poly2node = sparse(face2elem(:)*ones(1,3), face(:), 1, NP,N);
poly2node = (poly2node > 0);
NV = poly2node*ones(N, 1);
%centroid = poly2node*node./[NV, NV, NV];
%KP = pde.A(centroid);


% generate local IFE basis functions which are not associated with and DoFs
% IFEbasis(:,i,:,:) for minus subdomain (i=1) or plus domain (i=2)
iface = meshI.iface;
tface = meshI.tface;
iface2elem = meshI.iface2elem;
face2elemLoc = meshI.face2elemLoc;
PolyVolume = meshI.PolyVolume;
h = (PolyVolume(:,1) + PolyVolume(:,2)).^(1/3);

[iface2elemReduce, uniID] = unique(iface2elem);
% length(uniID) should be NP
ifaceReduce = iface(uniID,:);
[IntFnormal,IntFarea,unitNormal,unittgt1,unittgt2] = facenormaltgt(node,ifaceReduce);
[Fnormal,Farea,FunitNormal,Funittgt1,Funittgt2] = facenormaltgt(node,tface);
% IFEbasis contains the associated IFE function on each face
IFEbasis = zeros(size(face2elem,1),3,3);
% for two pieces of an ife function, only the normal vector is different
% on the minus subdomain (1) ife = n; on the plus subdomain (2) ife = n*bm/bp
IFEbasis(:,1,:) = unitNormal(face2elem,:);
IFEbasis(:,2,:) = unittgt1(face2elem,:);
IFEbasis(:,3,:) = unittgt2(face2elem,:);
piece1tmp = (face2elemLoc==1);
IFEbasis(~piece1tmp,1,:) = IFEbasis(~piece1tmp,1,:)*bm/bp;
BIFE =  zeros(length(face2elem),3);
% BIFE contains  (beta*grad vi.n)*|f|/3 for i=1,2,3
BIFE(:,1) = sum(squeeze(IFEbasis(:,1,:)).*Fnormal,2)/2;
BIFE(:,2) = sum(squeeze(IFEbasis(:,2,:)).*Fnormal,2)/2;
BIFE(:,3) = sum(squeeze(IFEbasis(:,3,:)).*Fnormal,2)/2;
BIFE(piece1tmp,:) = bm*BIFE(piece1tmp,:)/3.;
BIFE(~piece1tmp,:) = bp*BIFE(~piece1tmp,:)/3;
% generate the matrix for solving IFE projections which is a diagonal matrix
VdiagTmp1 = (bm*PolyVolume(:,1) + bm^2/bp*PolyVolume(:,2)).^(-1);
VdiagTmp2 = (bm*PolyVolume(:,1) + bp*PolyVolume(:,2)).^(-1);
VdiagTmp3 = (bm*PolyVolume(:,1) + bp*PolyVolume(:,2)).^(-1);
Vdiag1 = VdiagTmp1(face2elem);
Vdiag2 = VdiagTmp2(face2elem);
Vdiag3 = VdiagTmp3(face2elem);
Vdiag = zeros(3*length(face2elem),1);
Vdiag(1:3:end-2) = Vdiag1;
Vdiag(2:3:end-1) = Vdiag2;
Vdiag(3:3:end) = Vdiag3;
Mdiag = spdiags(Vdiag,0,3*length(face2elem),3*length(face2elem));
BIFE = reshape(Mdiag*reshape(BIFE',[],1),3,[])';
% this updated new BIFE(i,:) contains (old)BIFE(i,:)*M^(-1) where M is the
% coefficient matrix associated with the element for this face.
time = time + toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute projections to IFE spaces and form the matrices
idx2cube = meshI.idx2cube;
volume = meshI.tetVolume;
tetElem = meshI.tetElem;
tetElemLoc = meshI.tetElemLoc;

nnz = sum(NV.^2);
iiP = zeros(nnz, 1);
jjP = zeros(nnz, 1);
ssP = zeros(nnz, 1);
indexP = 0;
iiF = zeros(nnz, 1);
jjF = zeros(nnz, 1);
ssF = zeros(nnz, 1);
indexF = 0;
idx2cubeNew = -mesh.tLoc(idx2cube);
nnM = sum(NV(idx2cubeNew).^2);
iiM = zeros(nnM, 1);
jjM = zeros(nnM, 1);
ssM = zeros(nnM, 1);
indexM = 0;
unv = unique(NV);
hF = h(face2elem);
b = zeros(Ndof, 1);
ElemID = cell(1,2);
Bas1 = cell(1,2);
Bas2 = cell(1,2);
IFEC = cell(1,2);
ElemDoF = cell(1,2);
Xmpt = cell(1,2);

for kk = 1:length(unv) % group polys according to their # of nodes
    
    %% generate IFE functions including the IFE basis functions on each interface
    %% elements which are not associated with any DoFs, and projections of gradients
    %% and the IFE functions with matching face average
    tic
    nv = unv(kk);
    isCurrentFace = (NV(face2elem) == nv);
    CurrentFace = tface(isCurrentFace,:);
    CNF = size(CurrentFace,1);
    Currentfaceloc = face2elemLoc(isCurrentFace);
    Currentface2elem = face2elem(isCurrentFace);
    CurrentFnormal = FunitNormal(isCurrentFace,:);
    isCurrentElem = false(NP, 1);
    isCurrentElem(face2elem(isCurrentFace)) = true;
    CurrentPolyVolume = PolyVolume(isCurrentElem,:);
    % Current elem to node matrix
    currentPoly2node = poly2node(isCurrentElem, :);
    CNP = size(currentPoly2node, 1);
    currentPolyLocalIdx = zeros(NP, 1);
    currentPolyLocalIdx(isCurrentElem) = 1:CNP;
    
    % currentElem(i,:) contains the node index of the (current) i-th element
    [I, ~] = find(currentPoly2node');
    currentElem = reshape(I, nv, [])';
    % localIdx(i,j) contains the (current) i-th element having the
    % node(dof=1,2,3,4) on the j-th location (j is the global node index)
    localIdx = sparse(repmat((1:CNP)', 1, nv), currentElem, ones(CNP, 1)*(1:nv), CNP, N);
    clear currentPoly2node
    % Deal with current triangle face cases
    %isCTFace = isCurrentFace;
    tFace = face(isCurrentFace, :);
    NN = sum(isCurrentFace);
    subs1 = currentPolyLocalIdx(face2elem(isCurrentFace));
    subs2 = [ones(NN, 1); 2*ones(NN, 1); 3*ones(NN, 1)];
    val = BIFE(isCurrentFace,:);
    
    % compute the projection of gradients
    GP1 = zeros(CNP,3,nv);
    % coefficient of n t1 and t2 for the subelement on the subdomain 1
    for m = 1:3
        subs3 = full(localIdx((CurrentFace(:, m) - 1)*CNP + subs1));
        % subs1: new order (range in 1:CNP) of the element index
        % subs2: three components in the vector
        % subs3: the local index (DoF) (w.r.t. element) of the m-th node of this face
        GP1 = GP1 + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)], val(:), [CNP, 3, nv]);
    end
    GP2 = GP1; GP2(:,1,:) = bm/bp*GP2(:,1,:);
    % generate vector basis associated with interface plane
    unitNormalCurrent = unitNormal(isCurrentElem,:);
    unittgt1Current = unittgt1(isCurrentElem,:);
    unittgt2Current = unittgt2(isCurrentElem,:);
    NTCurrent = zeros(CNP,3,3);
    NTCurrent(:,:,1) = unitNormalCurrent;
    NTCurrent(:,:,2) = unittgt1Current;
    NTCurrent(:,:,3) = unittgt2Current;
    GP1 = mytimes(NTCurrent,GP1);
    GP2 = mytimes(NTCurrent,GP2);
    time = time + toc;
    % the updated GP contains three component values of the projection vector
    % generate IFE vectors on each face
    
    % generate some data structure on faces to compute IFE functions (not gradients)
    GP1tmp = zeros(NP,3,nv); GP2tmp = zeros(NP,3,nv);
    GP1tmp(isCurrentElem,:,:) = GP1; GP2tmp(isCurrentElem,:,:) = GP2;
    Gface = zeros(CNF,3,nv);
    Gface(Currentfaceloc==1,:,:) = GP1tmp(Currentface2elem(Currentfaceloc==1),:,:);
    Gface(Currentfaceloc==2,:,:) = GP2tmp(Currentface2elem(Currentfaceloc==2),:,:);
    clear GP1tmp GP2tmp
    % compute the IFE functions with the constant such that integral on
    % faces equal the virtual element function
    VPtmp1 = zeros(CNP,nv); % sum_i DoFi(u)|f|/3 on each face f
    VPtmp2 = zeros(CNP,nv); % grad*(xf-xm)|f|
    VPtmp3 = zeros(CNP,nv); % the surface area of each element repeated nv times
    FareaCurrent = Farea(isCurrentFace); %|f|/3
    for m = 1:3
        subs3 = full(localIdx((CurrentFace(:, m) - 1)*CNP + subs1));
        VPtmp1 = VPtmp1 + accumarray([subs1, subs3], FareaCurrent(:)/3, [CNP, nv]);
    end
    Xmf = (node(ifaceReduce(Currentface2elem,1),:) + node(ifaceReduce(Currentface2elem,2),:) +...
        node(ifaceReduce(Currentface2elem,3),:))/3; % the Xm point on the approximate interface plane
    Xf = (node(CurrentFace(:,1),:) + node(CurrentFace(:,2),:) +...
        node(CurrentFace(:,3),:))/3; % the centroid on each triangular face
    VPtmp2tmp = ((Xf - Xmf)).*FareaCurrent;
    for i = 1:nv
        tmp = sum(Gface(:,:,i).*VPtmp2tmp,2);
        VPtmp2(:,i) = accumarray(subs1, tmp,[CNP, 1]);
    end
    VPtmp3 = repmat(accumarray(subs1, FareaCurrent(:)),[1, nv]);
    IFEConst = (VPtmp1 - VPtmp2)./VPtmp3;
    clear currentPolyLocalIdx localIdx
    clear subs1 subs2 subs3
    
    %% form the stiffness matrices and stabilization matrices
    CurrentFace2DoF = zeros(NP,nv);
    CurrentFace2DoF(isCurrentElem,:) = currentElem;
    CurrentFace2DoF = CurrentFace2DoF(Currentface2elem,:);
    FunitNormalCurrent = FunitNormal(isCurrentFace,:);
    Funittgt1Current = Funittgt1(isCurrentFace,:);
    Funittgt2Current = Funittgt2(isCurrentFace,:);
    FNTCurrent = zeros(CNF,3,3);
    FNTCurrent(:,1,:) = Funittgt1Current;
    FNTCurrent(:,2,:) = Funittgt2Current;
    FNTCurrent(:,3,:) = FunitNormalCurrent;
    GfaceP = zeros(size(Gface));
    for n = 1:nv
        GfaceP(:,:,n) =  Gface(:,:,n) - sum(Gface(:,:,n).*CurrentFnormal,2).*CurrentFnormal;
    end
    GfaceP = mytimes(FNTCurrent,GfaceP);
    GfaceP = GfaceP(:,1:2,:);
    % generate vector basis associated with each face
    Fnode  = zeros(CNF,3,3);
    for l = 1:3
        Fnode(:,:,l) = node(CurrentFace(:,l),:);
    end
    FacePCoord = mytimes(FNTCurrent,Fnode);
    gradLambdatmp = zeros(CNF,2,3);  % barycentric coordinates
    gradLambda = zeros(CNF,2,nv);
    signtmp = ones(CNF,2); signtmp(:,1) = -1;
    gradLambdatmp(:,:,1) = (FacePCoord(:,[2,1],3) - FacePCoord(:,[2,1],2)).*signtmp./(2*Farea(isCurrentFace));
    gradLambdatmp(:,:,2) = (FacePCoord(:,[2,1],1) - FacePCoord(:,[2,1],3)).*signtmp./(2*Farea(isCurrentFace));
    gradLambdatmp(:,:,3) = -gradLambdatmp(:,:,1) - gradLambdatmp(:,:,2);
    for i = 1:nv
        for j = 1:3
            IDij = (CurrentFace(:,j) == CurrentFace2DoF(:,i));
            gradLambda(IDij,:,i) = gradLambdatmp(IDij,:,j);
        end
    end
    
    %     for n = 1:nv
    %         for m = 1:nv
    %             Xn = node(currentElem(:, n), :) - centroid(isCurrentElem,:);
    %             IminusP(:, n, m) = - 1/nv - dot(Xn, BP(:, :, m), 2)./volume(isCurrentElem); % stabilization
    %             if(n == m)
    %                 IminusP(:, n, m) = 1 + IminusP(:, n, m);
    %             end
    %         end
    %     end
    
    for n = 1:nv
        for m = 1:nv
            iiP(indexP+1:indexP + CNP) = currentElem(:, n);
            jjP(indexP+1:indexP + CNP) = currentElem(:, m);
            ssP(indexP+1:indexP + CNP) = bm*dot(GP1(:,:,n), GP1(:,:,m),2).*CurrentPolyVolume(:,1)+...
                bp*dot(GP2(:,:,n), GP2(:,:,m),2).*CurrentPolyVolume(:,2);
            indexP = indexP + CNP;
            
            if isfield(option,'scaling')
                scl = option.scl;
            else
                scl = 1;
            end
            
            iiF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, n);
            jjF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, m);
            ssF(indexF+1:indexF + CNF) = hF(isCurrentFace).^(scl).*...
                dot((gradLambda(:,:,m) - GfaceP(:,:,m)),(gradLambda(:,:,n) - GfaceP(:,:,n)),2).*...
                Farea(isCurrentFace);%.*hF(isCurrentFace).^(scl);
            indexF = indexF + CNF;
        end
    end
    
    %% compute the RHS with only evaluating the integral at the centroid
    %% approach 1 (will give the optimal order but slightly increase the error)
    %     currentElemSign = vSign(currentElem);
    %     currentElemSign1 = currentElemSign<=0; nv1 = sum(currentElemSign1,2);
    %     currentElemSign2 = currentElemSign>=0; nv2 = sum(currentElemSign2,2);
    %     currentElemNodeX = reshape(node(currentElem',1),size(currentElem,2),[])';
    %     currentElemNodeY = reshape(node(currentElem',2),size(currentElem,2),[])';
    %     currentElemNodeZ = reshape(node(currentElem',3),size(currentElem,2),[])';
    %     centroid1 = [sum(currentElemNodeX.*currentElemSign1,2)./nv1,...
    %         sum(currentElemNodeY.*currentElemSign1,2)./nv1, ...
    %         sum(currentElemNodeZ.*currentElemSign1,2)./nv1];
    %     centroid2 = [sum(currentElemNodeX.*currentElemSign2,2)./nv2,...
    %         sum(currentElemNodeY.*currentElemSign2,2)./nv2, ...
    %         sum(currentElemNodeZ.*currentElemSign2,2)./nv2];
    %     ft1 = CurrentPolyVolume(:,1).*pde.f(centroid1(:,1),centroid1(:,2),centroid1(:,3));
    %     ft2 = CurrentPolyVolume(:,2).*pde.f(centroid2(:,1),centroid2(:,2),centroid2(:,3));
    %
    %     Xm = (node(ifaceReduce(isCurrentElem,1),:) + node(ifaceReduce(isCurrentElem,2),:) +...
    %         node(ifaceReduce(isCurrentElem,3),:))/3;
    %     ifeEvalcen1 = zeros(CNP,nv); ifeEvalcen2 = zeros(CNP,nv);
    %     for n = 1:nv
    %         ifeEvalcen1(:,n) = sum(GP1(:,:,n).*(centroid1-Xm),2) + IFEConst(:,n);
    %         ifeEvalcen2(:,n) = sum(GP2(:,:,n).*(centroid2-Xm),2) + IFEConst(:,n);
    %     end
    %
    %     Currentb = ft1.*ifeEvalcen1 + ft2.*ifeEvalcen2;
    %     b = b + accumarray(currentElem(:), Currentb(:), [Ndof, 1]);
    %     clear ft isCurrentElem currentElem
    
    %% approach 2
    
    %idx2cubeNew = -mesh.tLoc(idx2cube);
    TetID = isCurrentElem(idx2cubeNew);
    currentvolume = volume(TetID);
    X1 = node(tetElem(TetID,1),:);
    X2 = node(tetElem(TetID,2),:);
    X3 = node(tetElem(TetID,3),:);
    X4 = node(tetElem(TetID,4),:);
    ng = 1;
    [gx, gy, gz] = gaussPtetra(X1, X2, X3, X4, ng);
    gw = gaussWtetra(ng);
    Xm = (node(ifaceReduce(isCurrentElem,1),:) + node(ifaceReduce(isCurrentElem,2),:) +...
        node(ifaceReduce(isCurrentElem,3),:))/3;
    
    Bastmp1 = zeros(size(isCurrentElem,1),3,nv); Bastmp2 = Bastmp1;
    Bastmp1(isCurrentElem,:,:) = GP1;
    Bastmp2(isCurrentElem,:,:) = GP2;
    Bastmp1 = Bastmp1(idx2cubeNew(TetID),:,:);
    Bastmp2 = Bastmp2(idx2cubeNew(TetID),:,:);
    Bas = zeros(sum(TetID),3,nv);
    %%%
    IFEc0 = zeros(size(isCurrentElem,1),nv);
    IFEc0(isCurrentElem,:) = IFEConst;
    IFEc0 = IFEc0(idx2cubeNew(TetID),:);
    %%%
    Xmtet = zeros(size(isCurrentElem,1),3);
    Xmtet(isCurrentElem,:) = Xm;
    Xmtet = Xmtet(idx2cubeNew(TetID),:);
    %%%
    currentTetDoF = zeros(size(isCurrentElem,1),nv);
    currentTetDoF(isCurrentElem,:) = currentElem;
    currentTetDoF = currentTetDoF(idx2cubeNew(TetID),:);
    %%%
    piecetmp = tetElemLoc(TetID);
    Bas(piecetmp==1,:,:) = Bastmp1(piecetmp==1,:,:);
    Bas(piecetmp==2,:,:) = Bastmp2(piecetmp==2,:,:);
    if option.rhs == 1
        ft = pde.f(gx,gy,gz);
        Currentb = zeros(size(ft));
        cnt = size(gx,1);
        bii = zeros(nv*cnt,1);
        count = 0 ;
        for i = 1:nv
            Basi = squeeze(Bas(:,:,i)); IFEc0i = IFEc0(:,i);
            uhi = (gx-Xmtet(:,1)).*Basi(:,1) + (gy-Xmtet(:,2)).*Basi(:,2) +...
                (gz-Xmtet(:,3)).*Basi(:,3) + IFEc0i;
            Currentb(:,i) = sum(gw*ft.*uhi,2).*currentvolume;
            bii(count+1:count+cnt) = currentTetDoF(:,i);
            count = count+cnt;
        end
        b = b + sparse(bii,1,reshape(Currentb,[],1),Ndof,1);
    end
    
    if option.mass == 1 % generate mass matrix
        cnt = size(gx,1);
        for i = 1:nv
            for j = 1:nv
                iiM(indexM+1:indexM + cnt) = currentTetDoF(:, i);
                jjM(indexM+1:indexM + cnt) = currentTetDoF(:, j);
                
                Basi = squeeze(Bas(:,:,i)); IFEc0i = IFEc0(:,i);
                Basj = squeeze(Bas(:,:,j)); IFEc0j = IFEc0(:,j);
                uhi = (gx-Xmtet(:,1)).*Basi(:,1) + (gy-Xmtet(:,2)).*Basi(:,2) +...
                    (gz-Xmtet(:,3)).*Basi(:,3) + IFEc0i;
                uhj = (gx-Xmtet(:,1)).*Basj(:,1) + (gy-Xmtet(:,2)).*Basj(:,2) +...
                    (gz-Xmtet(:,3)).*Basj(:,3) + IFEc0j;
                ssM(indexM+1:indexM + cnt) = sum(gw*uhi.*uhj,2).*currentvolume;
                indexM = indexM+cnt;
            end
        end
        
    end

ElemID{kk} = isCurrentElem;
Bas1{kk} = GP1; Bas2{kk} = GP2;
IFEC{kk} = IFEConst;
ElemDoF{kk} = currentElem;
Xmpt{kk} = Xm;
end
K = sparse(iiP, jjP, ssP, Ndof, Ndof);
S = sparse(iiF, jjF, ssF, Ndof, Ndof);

if option.mass == 1
    M = sparse(iiM, jjM, ssM, Ndof, Ndof);
    femI.M = M;
end

femI.Bstime = time;
femI.Bas1 = Bas1; femI.Bas2 = Bas2;
femI.IFEconst = IFEC;
femI.K = K; femI.S = S; 
femI.b = b;
femI.ElemID = ElemID;
femI.ElemDoF = ElemDoF;
femI.Xmpt = Xmpt;
femI.NV = unv;
