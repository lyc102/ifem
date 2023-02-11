%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. use vectorization to generate face and face2elem data structure
%% this vectorization code does not work for uniform tets

tic

Np = size(mesh.p,1);
node = [mesh.p;mesh.eIntP];
allCutElem = mesh.t(mesh.tLoc<0,:);
isCutElem = (mesh.tLoc<0);
vSign = [mesh.pLoc;zeros(size(mesh.eIntP,1),1)];
isInterfaceNode = false(size(node,1), 1);
isInterfaceNode(allCutElem(:)) = true; % include the vertices of allCutElem 
numCutElemV1 = sum(isInterfaceNode);
isInterfaceNode(Np+1:end) = true;    % and the cut points and the aux points
allCutElemReoderTmp = zeros(size(node,1),1);
numCutElemV2 = sum(isInterfaceNode);
allCutElemReoderTmp(isInterfaceNode) = 1:numCutElemV2;
allCutElemReoder = allCutElemReoderTmp(allCutElem);
interfaceNode = node(isInterfaceNode,:);
ntI = -min(mesh.tLoc);
intID = find(mesh.tLoc<0);

tetElem = zeros(12*ntI,4);
tetElemLoc = zeros(12*ntI,1);
idx2cube = zeros(12*ntI,1);
tetcount = 0;

for i = 1:ntI
    tID = intID(i);
    t_e = mesh.t_e(tID,:);
    pLocK = mesh.pLoc(mesh.t(tID,:));
    vert0 = mesh.p(mesh.t(tID,:),:);
    nodeid = [allCutElemReoder(i,:),numCutElemV1-mesh.eLoc(t_e(mesh.eLoc(t_e)<0))'];
    vert2 = vert0(pLocK>0,:); % plus domain
    vert1 = vert0(pLocK<0,:); % minus domain
    intpt0 = mesh.eIntP(-mesh.eLoc(t_e(mesh.eLoc(t_e)<0)),:);
    id1 = [find(pLocK<0)', 5:4+size(intpt0,1)];
    id2 = [find(pLocK>0)', 5:4+size(intpt0,1)];
    p1 = [vert1;intpt0]; DT = delaunayTriangulation(p1); t1 = DT.ConnectivityList;
    p2 = [vert2;intpt0]; DT = delaunayTriangulation(p2); t2 = DT.ConnectivityList;
    tetElem(tetcount+1:tetcount+size(t1,1),:) = nodeid(id1(t1));
    idx2cube(tetcount+1:tetcount+size(t1,1)) = tID;
    tetElemLoc(tetcount+1:tetcount+size(t1,1)) = 1; % inside subdomain
    tetcount = tetcount+size(t1,1);
    tetElem(tetcount+1:tetcount+size(t2,1),:) = nodeid(id2(t2));
    idx2cube(tetcount+1:tetcount+size(t2,1)) = tID;
    tetElemLoc(tetcount+1:tetcount+size(t2,1)) = 2; % outside subdomain
    tetcount = tetcount+size(t2,1);
end
tetElem(tetcount+1:end,:) = [];
tetElemLoc(tetcount+1:end,:) = [];
idx2cube(tetcount+1:end,:) = [];
[tetElem, ortidx, volume] = fixorder3(interfaceNode, tetElem);
% there are some tets which are coplane, need to get rid of them
idPlane = (volume<=10^(-12));
tetElem(idPlane,:) = [];
volume(idPlane,:) = [];
tetElemLoc(idPlane,:) = [];
idx2cube(idPlane,:) = [];
localidx2globalidx = find(isInterfaceNode); % tetElem points to interfaceNode
tetElem = localidx2globalidx(tetElem);  
% reorder tetElem such that the index of its nodes are local, i.e.,
% starting from 1,2,...

% there are some tets which are on the interface, but contained in the surrounding tets
% need to get rid of them (but the following algorithm may not be robust)
tetSign = vSign(tetElem);
idSliver = (sum(abs(tetSign),2)==0);
tetElem(idSliver,:) = [];
volume(idSliver,:) = [];
tetElemLoc(idSliver,:) = [];
idx2cube(idSliver,:) = [];
PolyVolume = accumarray([idx2cube,tetElemLoc],volume);
PolyVolume = PolyVolume(mesh.tLoc<0,:); 
h = (PolyVolume(:,1) + PolyVolume(:,2)).^(1/3);
% the first component contains the volume of the polyhedron on the
% subdomain 1, the first component contains the volume of the polyhedron
% on the subdomain 2

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get triangular faces for interrior elements and interface

tic
localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
NT = size(tetElem,1);
tface = zeros(4*NT, 3);
tface2elem = zeros(4*NT, 1);
iface = zeros(4*NT, 3);
iface2elem = zeros(4*NT, 1);
% find the interior tet elements
isIntTet = min(vSign(tetElem),[], 2) == -1; % can not be == -1 as there is sliver
intTet = tetElem(isIntTet,:);
% find the corresponding cube indices
intIdx2cube = idx2cube(isIntTet); 
% find triangular interface
T = auxstructure3(intTet);
neighbor = T.neighbor; % if a face is on the boundary, then the neighbor element index is itself
clear T;
tmp = (1:size(intTet, 1))';
ct = 0;
ci = 0;
for i = 1:4
    face = intTet(:, localFace(i,:));
    % find the triangle faces of polyhedron
    % 1. face and its neighbor associated to different cubes, and
    % 2. face is not on the boundary of all cut elements (which are squares)
    isPolyTriFace = ((neighbor(:, i) == tmp) | (intIdx2cube ~= intIdx2cube(neighbor(:,i)))) &...
        (sum(abs(vSign(face)), 2) > 0);% & (sum(abs(vSign(face)), 2) > 0);
    c2 = ct + sum(isPolyTriFace);
    tface((ct+1):c2,:) = face(isPolyTriFace,:);
    tface2elem((ct+1):c2,:) = intIdx2cube(isPolyTriFace); % the indices of the polyhedron
    ct = c2;
    % note that only interior elements are treated
    
    % find the triangle faces on interface with normal points to exterior
    % 1. face is on the boundary and
    % 2. all vertices are on the interface
    isInterface = (sum(abs(vSign(face)), 2) == 0);
    
    % add to interface    
    c4 = ci + sum(isInterface);
    iface((ci+1):c4,:) = face(isInterface,:); % the interface tri faces.
    iface2elem((ci+1):c4,:) = intIdx2cube(isInterface); % the indices of the polyhedron
    ci = c4;
end
iface((ci+1):end,:) = [];
iface2elem((ci+1):end,:) = [];
c2old = c2;
% plot to check
%trisurf(tface(1:c2,:),node(:,1),node(:,2),node(:,3))
%trisurf(iface(1:ci,:),node(:,1),node(:,2),node(:,3))

% Find the triangular faces for exterior elements 
extTet = tetElem(~isIntTet, :);
extIdx2cube = idx2cube(~isIntTet);

T = auxstructure3(extTet);
neighbor = T.neighbor;
clear T;
tmp = (1:size(extTet,1))';
for i = 1:4
    face = extTet(:, localFace(i,:));
    % find the triangle faces of polyhedron
    % 1. face and its neighbor associated to different cubes, and
    % 2. face is not on the boundary of all cut elements (which are squares)
    isPolyTriFace = ((neighbor(:, i) == tmp) | (extIdx2cube ~= extIdx2cube(neighbor(:,i)))) &...
        (sum(abs(vSign(face)), 2) > 0);
    c2 = ct + sum(isPolyTriFace);
    tface((ct+1):c2,:) = face(isPolyTriFace,:);
    tface2elem((ct+1):c2,:) = extIdx2cube(isPolyTriFace); 
    % index of exterior polyhedron is append to the end of elem
    ct = c2;
end
toc

tface((ct+1):end,:) = [];    % remove empty meomory
tface2elem((ct+1):end) = [];
face2elemLoc = [ones(c2old,1);2*ones(c2-c2old,1)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate data structure for computing projection of Hdiv functions
face2elem = tface2elem;
N = size(node, 1);  
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

[NewFace,NewFaceID1,NewFaceID2] = unique(sort(tface,2),'rows');
% NewFaceID2 containing some index with multiple times which should be 
% the DoF of the faces in tface
% NewFace is already sorted
% this means the order of the all newly added faces including those
% non-interfaces should follow the sorted one NewFace
NewFaceNintID = find(abs(sum(vSign(NewFace),2))==3);
NewFaceNint = NewFace(NewFaceNintID,:);
OldFaceID = unique(reshape(mesh.t_f(mesh.tLoc<0,:),[],1));
OldFace = mesh.f(OldFaceID,:);
OldFtmp = abs(sum(vSign(OldFace),2))==3;
OldFaceNintID = OldFaceID(OldFtmp); 
OldFaceNint = OldFace(OldFtmp,:);
[OldFaceNintSort,OldFaceNintSortID] = sortrows(sort(OldFaceNint,2));
% dF = NewFaceNint - OldFaceNintSort; this should be zero

BasicFace  = true(size(mesh.f,1),1); BasicFace(OldFaceID) = false;
BasicFaceNum = sum(BasicFace);
TotalOldFace = zeros(size(mesh.f,1),1);
TotalOldFace(BasicFace) = 1:BasicFaceNum;
IDtmp = OldFaceNintID(OldFaceNintSortID);
TotalOldFace(IDtmp) = BasicFaceNum + NewFaceNintID;

NFdof = size(mesh.f,1) + 2*sum(mesh.fLoc<0);
NumFace = NFdof;
Fgdof = zeros(NumFace,3);  % contains index of nodes
Fgdof(1:BasicFaceNum,:) = mesh.f(BasicFace,:);
Fgdof(BasicFaceNum+1:end,:) = NewFace;
NintElem = mesh.tLoc>0;
g2lFdofNint = TotalOldFace(mesh.t_f(NintElem,:));
face2FDoF = BasicFaceNum+NewFaceID2;

poly2face = sparse(face2elem, face2FDoF, 1, NP,NFdof);
poly2face = (poly2face > 0);
NF = poly2face*ones(NFdof, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate local IFE basis functions which are not associated with any DoFs
% IFEbasis(:,i,:,:) for minus subdomain (i=1) or plus domain (i=2)
[iface2elemReduce, uniID] = unique(iface2elem);
% length(uniID) should be NP
ifaceReduce = iface(uniID,:);
[IntFnormal,IntFarea,unitNormal,unittgt1,unittgt2] = facenormaltgt(node,ifaceReduce);
[Fnormal,Farea,FunitNormal,Funittgt1,Funittgt2] = facenormaltgt(node,tface);
vOriet1 = [0,0,1]; vOriet2 = [1,0,0]; vOriet3 = [0,1,0];
vid1 = sum(FunitNormal.*vOriet1,2)>10^(-12);
vid2 = (abs(sum(FunitNormal.*vOriet1,2))<=10^(-12) & sum(FunitNormal.*vOriet2,2)>10^(-12));
vid3 = (abs(sum(FunitNormal.*vOriet1,2))<=10^(-12) & abs(sum(FunitNormal.*vOriet2,2))<10^(-12)...
    & sum(FunitNormal.*vOriet3,2)>=10^(-12));
vid = (vid1+vid2+vid3>0);
face2FDoFSign = -ones(size(tface,1),1);
face2FDoFSign(vid) = 1;
% IFEbasis contains the associated IFE function on each face
% IFEbasisA is for curl term and IFEbasisB is for u term
IFEbasisA = zeros(size(face2elem,1),3,3);
% for two pieces of an ife function, only the normal vector is different
% on the minus subdomain (1) ife = n; on the plus subdomain (2) ife = n*bm/bp
IFEbasisA(:,1,:) = unitNormal(face2elem,:);
IFEbasisA(:,2,:) = unittgt1(face2elem,:);
IFEbasisA(:,3,:) = unittgt2(face2elem,:);
piece1tmp = (face2elemLoc==1);
IFEbasisA(~piece1tmp,2,:) = IFEbasisA(~piece1tmp,2,:)*am/ap;
IFEbasisA(~piece1tmp,3,:) = IFEbasisA(~piece1tmp,3,:)*am/ap;
% Xmf the centroid of a triangular face on the interface of each interface
% element but assigned to each face
% Xf the centroid of each triangular face
Xmf = (node(ifaceReduce(face2elem,1),:) + node(ifaceReduce(face2elem,2),:) +...
    node(ifaceReduce(face2elem,3),:))/3; % the Xm point on the approximate interface plane
Xf = (node(tface(:,1),:) + node(tface(:,2),:) +...
    node(tface(:,3),:))/3;

%% generate boundary integral for computing projections
% for the curl term, need to compute surface curl (rot) on each face
BIFEA =  zeros(length(face2elem),3);
% BIFEA contains  (alpha*(grad vi).(xf -xm))for i=1,2,3 (vi is the test function)
BIFEA(:,1) = sum(squeeze(IFEbasisA(:,1,:)).*(Xf - Xmf),2); 
BIFEA(:,2) = sum(squeeze(IFEbasisA(:,2,:)).*(Xf - Xmf),2); 
BIFEA(:,3) = sum(squeeze(IFEbasisA(:,3,:)).*(Xf - Xmf),2); 
% no need to mutiply |F| as the DoF is integral
BIFEA(piece1tmp,:) = am*BIFEA(piece1tmp,:);
BIFEA(~piece1tmp,:) = ap*BIFEA(~piece1tmp,:);
% BIFEB contains  (beta*(curl vi)times(xf -xm))times n for i=1,2,3
% position to match the structure of mytimes
%%%%
% mutiply the inverse the matrix on the left hand side of each interface
% element for computing projection
% this is a blocal diagonal matrix
% for BIFEA (curl):
VdiagTmp1 = (am*PolyVolume(:,1) + ap*PolyVolume(:,2)).^(-1);
VdiagTmp2 = (am*PolyVolume(:,1) + am^2/ap*PolyVolume(:,2)).^(-1);
VdiagTmp3 = (am*PolyVolume(:,1) + am^2/ap*PolyVolume(:,2)).^(-1);
Vdiag1 = VdiagTmp1(face2elem);
Vdiag2 = VdiagTmp2(face2elem);
Vdiag3 = VdiagTmp3(face2elem);
Mdiag = zeros(size(Vdiag1,1),3,3);
Mdiag(:,1,1) = Vdiag1; Mdiag(:,2,2) = Vdiag2; Mdiag(:,3,3) = Vdiag3;
BIFEA = mytimes(Mdiag,BIFEA);
clear IFEbasisA FNT VdiagTmp1 VdiagTmp2 VdiagTmp3 Vdiag1 Vdiag2 Vdiag3 Mdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate IFE function selves
IFEbasisA = zeros(size(tetElem,1),3,3);
idx2cubeNew = -mesh.tLoc(idx2cube);
IFEbasisA(:,1,:) = unitNormal(idx2cubeNew,:);
IFEbasisA(:,2,:) = unittgt1(idx2cubeNew,:);
IFEbasisA(:,3,:) = unittgt2(idx2cubeNew,:);
piece1tmp = (tetElemLoc==1);
IFEbasisA(~piece1tmp,2,:) = IFEbasisA(~piece1tmp,2,:)*am/ap;
IFEbasisA(~piece1tmp,3,:) = IFEbasisA(~piece1tmp,3,:)*am/ap;

TetNode1 = node(tetElem(:,1),:); TetNode2 = node(tetElem(:,2),:);
TetNode3 = node(tetElem(:,3),:); TetNode4 = node(tetElem(:,4),:);
Xm = (node(ifaceReduce(:,1),:) + node(ifaceReduce(:,2),:) +...
        node(ifaceReduce(:,3),:))/3;
Xmt = Xm(idx2cubeNew,:);
centroid = (TetNode1 + TetNode2 + TetNode3 + TetNode4)/4;
BIFEPhi =  zeros(size(tetElem,1),3);
ElemVolume = PolyVolume(:,1) + PolyVolume(:,2);
BIFEPhi(:,1) = sum(squeeze(IFEbasisA(:,1,:)).*(centroid - Xmt),2).*...
    volume./ElemVolume(idx2cubeNew); 
BIFEPhi(:,2) = sum(squeeze(IFEbasisA(:,2,:)).*(centroid - Xmt),2).*...
    volume./ElemVolume(idx2cubeNew); 
BIFEPhi(:,3) = sum(squeeze(IFEbasisA(:,3,:)).*(centroid - Xmt),2).*...
    volume./ElemVolume(idx2cubeNew); 
BIFEPhi(piece1tmp,:) = am*BIFEPhi(piece1tmp,:);
BIFEPhi(~piece1tmp,:) = ap*BIFEPhi(~piece1tmp,:);
VdiagTmp1 = (am*PolyVolume(:,1) + ap*PolyVolume(:,2)).^(-1);
VdiagTmp2 = (am*PolyVolume(:,1) + am^2/ap*PolyVolume(:,2)).^(-1);
VdiagTmp3 = (am*PolyVolume(:,1) + am^2/ap*PolyVolume(:,2)).^(-1);
Vdiag1 = VdiagTmp1(idx2cubeNew);
Vdiag2 = VdiagTmp2(idx2cubeNew);
Vdiag3 = VdiagTmp3(idx2cubeNew);
Mdiag = zeros(size(Vdiag1,1),3,3);
Mdiag(:,1,1) = Vdiag1; Mdiag(:,2,2) = Vdiag2; Mdiag(:,3,3) = Vdiag3;
BIFEPhi = mytimes(Mdiag,BIFEPhi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFNT = zeros(length(ifaceReduce),3,3);
IFNT(:,:,1) = unitNormal;
IFNT(:,:,2) = unittgt1;
IFNT(:,:,3) = unittgt2;

nnz = sum(NF.^2);
iiP = zeros(nnz, 1);
jjP = zeros(nnz, 1);
ssKC = zeros(nnz, 1);
indexP = 0;
iiF = zeros(nnz, 1);
jjF = zeros(nnz, 1);
ssSC = zeros(nnz, 1);
indexF = 0;
unf = unique(NF);
hF = h(face2elem);
b = zeros(NFdof, 1);

for kk = 1:length(unf) % group polys according to their # of nodes
    
    tic
    %% generate IFE functions including the IFE basis functions on each interface
    %% elements which are not associated with any DoFs, and projections of gradients
    %% and the IFE functions with matching face average
    nf = unf(kk);
    % some data related to element
    isCurrentFace = (NF(face2elem) == nf);
    CurrentFace = tface(isCurrentFace,:);
    CurrentFaceF = face2FDoF(isCurrentFace,:);
    CurrentFaceFS = face2FDoFSign(isCurrentFace,:);
    CNF = size(CurrentFaceF,1);
    Currentfaceloc = face2elemLoc(isCurrentFace);
    Currentface2elem = face2elem(isCurrentFace);
    CurrentFnormal = FunitNormal(isCurrentFace,:);
    FareaCurrent = Farea(isCurrentFace);
    % some data related to face
    isCurrentElem = false(NP, 1);
    isCurrentElem(face2elem(isCurrentFace)) = true;
    CurrentPolyVolume = PolyVolume(isCurrentElem,:);
    IFNTCurrent = IFNT(isCurrentElem,:,:);
    isCurrentTet = isCurrentElem(idx2cubeNew);
    % Current elem to node matrix
    currentPoly2face = poly2face(isCurrentElem, :);
    CNP = size(currentPoly2face, 1);
    currentPolyLocalIdx = zeros(NP, 1);
    currentPolyLocalIdx(isCurrentElem) = 1:CNP;
    
    % currentElem(i,:) contains the edge index of the (current) i-th element
    [I, ~] = find(currentPoly2face');
    currentElem = reshape(I, nf, [])';
    % localIdx(i,j) contains the (current) i-th element having the
    % node(dof=1,2,3,4) on the j-th location (j is the global node index)
    localIdx = sparse(repmat((1:CNP)', 1, nf), currentElem, ones(CNP, 1)*(1:nf), CNP, NFdof);
    clear currentPoly2face
    
    %% Deal with current triangle face cases
    %% (1) for projection of u in Hdiv
    subs1 = currentPolyLocalIdx(face2elem(isCurrentFace));
    subs2 = [ones(CNF, 1); 2*ones(CNF, 1); 3*ones(CNF, 1)];
    val = BIFEA(isCurrentFace,:); 
    % compute the projection of gradients
    %GPA1 = zeros(CNP,3,nf);
    % coefficient of n t1 and t2 for the subelement on the subdomain 1
    subs3 = full(localIdx((CurrentFaceF(:, 1) - 1)*CNP + subs1));
    SignCurrentF = repmat(CurrentFaceFS,3,1);
    % subs1: new order (range in 1:CNP) of the element index
    % subs2: three components in the vector
    % subs3: the local index (DoF) (w.r.t. element) of this face
    GPA1 = accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)], SignCurrentF.*val(:), [CNP, 3, nf]);
    % the updated GP contains three component values of the projection vector   
    % generate IFE vectors on each face
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute int_K div(vh) varphi_h dx
    % compute div(vh) first
    divu = zeros(NP,nf);
    divu(isCurrentElem,:) = accumarray([subs1,subs3], CurrentFaceFS, [CNP, nf]);
    divuTet = divu(idx2cubeNew(isCurrentTet),:);
    val = BIFEPhi(idx2cubeNew(isCurrentTet),:);
    [unitmp,~,idtmp] = unique(idx2cubeNew(isCurrentTet));
    subs4 = 1:length(unitmp); subs4 = subs4(idtmp);
    for m = 1:nf
        for i = 1:3
        GPA1(:,i,m) = GPA1(:,i,m) - accumarray(subs4',divuTet(:,m).*val(:,i));
        end
    end
    
    GPA2 = GPA1; GPA2(:,2,:) = am/ap*GPA2(:,2,:); GPA2(:,3,:) = am/ap*GPA2(:,3,:);
    GPA1 = mytimes(IFNTCurrent,GPA1);
    GPA2 = mytimes(IFNTCurrent,GPA2);
    
    clear currentPolyLocalIdx localIdx
    clear subs1 subs2 subs3
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CurrentFace2DoF = zeros(NP,nf);
    CurrentFace2DoF(isCurrentElem,:) = currentElem;
    CurrentFace2DoF = CurrentFace2DoF(Currentface2elem,:);
    % CurrentFace2DoF contains the DoFs of the element associated with each face
    
    % compute curl uh.n on each face (uh is an IFE function)
    GPA1tmp = zeros(NP,3,nf); GPA2tmp = zeros(NP,3,nf);
    GPA1tmp(isCurrentElem,:,:) = GPA1; GPA2tmp(isCurrentElem,:,:) = GPA2;
    GAface = zeros(CNF,3,nf);
    GAface(Currentfaceloc==1,:,:) = GPA1tmp(Currentface2elem(Currentfaceloc==1),:,:);
    GAface(Currentfaceloc==2,:,:) = GPA2tmp(Currentface2elem(Currentfaceloc==2),:,:);
    clear GPA1tmp GPA2tmp     
    GAfaceP = zeros(size(GAface,1),nf);
    for n = 1:nf
        GAfaceP(:,n) =  sum(squeeze(GAface(:,:,n)).*CurrentFnormal,2);
    end
    
    % compute u.n uh on each face by the DoFs
    TrueUNTmp = CurrentFaceFS;
    TrueUN = zeros(CNF,nf);
    for i = 1:nf
        IDi = (CurrentFaceF == CurrentFace2DoF(:,i));
        TrueUN(IDi,i) = TrueUNTmp(IDi);
    end
    
    for n = 1:nf
        for m = 1:nf
            iiP(indexP+1:indexP + CNP) = currentElem(:, n);
            jjP(indexP+1:indexP + CNP) = currentElem(:, m);
            ssKU(indexP+1:indexP + CNP) = am*dot(GPA1(:,:,n), GPA1(:,:,m),2).*CurrentPolyVolume(:,1)+...
                ap*dot(GPA2(:,:,n), GPA2(:,:,m),2).*CurrentPolyVolume(:,2);
            indexP = indexP + CNP;
            
            iiF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, n);
            jjF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, m);
            ssSU(indexF+1:indexF + CNF) = dot(TrueUN(:,m) - GAfaceP(:,m).*Farea(isCurrentFace),...
                TrueUN(:,n) - GAfaceP(:,n).*Farea(isCurrentFace),2).*Farea(isCurrentFace).^(-1).*hF(isCurrentFace);
            indexF = indexF + CNF;
        end
    end
    
end

KU = sparse(iiP, jjP, ssKU, NFdof, NFdof);
SU = sparse(iiF, jjF, ssSU, NFdof, NFdof);

MatUI = KU + SU;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
localFace = [1 3 2; 1 2 4; 1 4 3; 2,3,4];
elem2faceSign = zeros(sum(mesh.tLoc>0),4);
for i = 1:4
    v1 = node(mesh.t(mesh.tLoc>0,localFace(i,2)),:) - node(mesh.t(mesh.tLoc>0,localFace(i,1)),:);
    v2 = node(mesh.t(mesh.tLoc>0,localFace(i,3)),:) - node(mesh.t(mesh.tLoc>0,localFace(i,1)),:);
    OneFaceNormal = cross(v1,v2);
    vOriet1 = [0,0,1]; vOriet2 = [1,0,0];
    vid1 = sum(OneFaceNormal.*vOriet1,2)>10^(-12);
    vid2 = (abs(sum(OneFaceNormal.*vOriet1,2))<=10^(-12) & sum(OneFaceNormal.*vOriet2,2)>10^(-12));
    vid3 = (abs(sum(OneFaceNormal.*vOriet1,2))<=10^(-12) & abs(sum(OneFaceNormal.*vOriet2,2))<10^(-12)...
        & sum(OneFaceNormal.*vOriet3,2)>=10^(-12));
    vid = (vid1+vid2+vid3>0);
    OneFaceNormalSign = -ones(sum(mesh.tLoc>0),1);
    OneFaceNormalSign(vid) = 1;
    elem2faceSign(:,i) = OneFaceNormalSign;
end

[Dlambda,volume,elemSign] = gradbasis3(node,mesh.t(mesh.tLoc>0,:));
elem2face = g2lFdofNint;%==(:,[4,3,2,1]);


tuv = [1,1,1]; TUV = repmat(tuv,size(elem2face,1),1);
%localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2]; % ascend ordering
MatNI = sparse(NFdof,NFdof);
rhsdiv = sparse(NFdof,1);
for i = 1:4
    for j = i:4
        % local to global index map
        ii = double(elem2face(:,i));
        jj = double(elem2face(:,j));
        i1 = localFace(i,1); i2 = localFace(i,2); i3 = localFace(i,3);
        j1 = localFace(j,1); j2 = localFace(j,2); j3 = localFace(j,3);
        % computation of mass matrix --- (phi_i, phi_j)
        Mij = 1/20*volume*4.*( ...
            (1+(i1==j1))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
            mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
            +(1+(i1==j2))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
            mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
            +(1+(i1==j3))*dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
            mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
            +(1+(i2==j1))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
            mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
            +(1+(i2==j2))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
            mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
            +(1+(i2==j3))*dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
            mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)...
            +(1+(i3==j1))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
            mycross(Dlambda(:,:,j2),Dlambda(:,:,j3),2),2)...
            +(1+(i3==j2))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
            mycross(Dlambda(:,:,j3),Dlambda(:,:,j1),2),2)...
            +(1+(i3==j3))*dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
            mycross(Dlambda(:,:,j1),Dlambda(:,:,j2),2),2)).*...
            elem2faceSign(:,i).*elem2faceSign(:,j);
            %elem2faceArea(:,i).^(-1).*elem2faceArea(:,j).^(-1).*...
        if (j==i)
            MatNI = MatNI + sparse(ii,jj,Mij,NFdof,NFdof);
        else
            MatNI = MatNI + sparse([ii;jj],[jj;ii],[Mij; Mij],NFdof,NFdof);
        end
    end
%     ri = 1/20*volume*2.*( ...
%         dot(mycross(Dlambda(:,:,i2),Dlambda(:,:,i3),2), ...
%         TUV,2)...
%         +dot(mycross(Dlambda(:,:,i3),Dlambda(:,:,i1),2), ...
%         TUV,2)...
%         +dot(mycross(Dlambda(:,:,i1),Dlambda(:,:,i2),2), ...
%         TUV,2)).*elem2faceSign(:,i);
%     rhsdiv = rhsdiv + sparse(ii,1,ri,NFdof,1);
end

%rhsdiv = full(rhsdiv);

MassTotal = MatNI + MatUI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% a test for solving the linear system

b = ones(size(MassTotal,1),1);
% alpha = 10;
% L = ichol(MassTotal,struct('type','ict','droptol',1e-8,'diagcomp',alpha));
% [u,flag,relres,iter,resvec] = pcg(MassTotal,b,1e-8,2*10^4,L,L');

% D=diag(diag(MassTotal).^(-1/2));
% alpha = 10;
% L = ichol(D*MassTotal*D,struct('type','ict','droptol',1e-8,'diagcomp',alpha));
% [u,flag,relres,iter,resvec] = pcg(D*MassTotal*D,D*b,1e-8,2*10^4,L,L');

D = diag(MassTotal(1:BasicFaceNum,1:BasicFaceNum));
M = MassTotal(BasicFaceNum+1:end,BasicFaceNum+1:end);
x0 = zeros(size(b));
[x,info] = blockpcg(MassTotal,b,x0,1e-8,10^3,D,M);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test that the projection matrix will preseve a constant vector
% [~,Allarea,AllunitNormal,~,~] = facenormaltgt(node,Fgdof);
% vOriet1 = [0,0,1]; vOriet2 = [1,0,0];
% vid1 = sum(AllunitNormal.*vOriet1,2)>10^(-12);
% vid2 = (abs(sum(AllunitNormal.*vOriet1,2))<=10^(-12) & sum(AllunitNormal.*vOriet2,2)>10^(-12));
% vid3 = (abs(sum(AllunitNormal.*vOriet1,2))<=10^(-12) & abs(sum(AllunitNormal.*vOriet2,2))<10^(-12)...
%     & sum(AllunitNormal.*vOriet3,2)>=10^(-12));
% vid = (vid1+vid2+vid3>0);
% Allface2FDoFSign = -ones(size(Fgdof,1),1);
% Allface2FDoFSign(vid) = 1;
% tu = Allarea.*(AllunitNormal*tuv').*Allface2FDoFSign;
% 
% for ElemID = 1:size(GPA1,1)
%     tuv_test = squeeze(GPA1(ElemID,:,:))*tu(currentElem(ElemID,:));
%     if norm(tuv_test-tuv')>10^(-10)
%         stp=1;
%     end
% end




%[bc,mapper] = boundaryFace3D(node,Fgdof,[1,1,1,1,1,1]);
%A = MatU(mapper,mapper);