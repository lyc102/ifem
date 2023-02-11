%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

domain = [-1,1,-1,1,-1,1];
bc = [1,1,1,1,1,1];
nx = 10;
ny = nx;
nz = nx;
r1 = 1; r2 = 0; r3 = 0; x0 = pi/10^4; y0 = 0; z0 = 0;
bm = 1; bp = 1; am = 1; ap = 1;
pde = ConstFun(am,ap,bm,bp,r1,r2,r3,x0,y0,z0);
mesh = genMesh3D(domain, nx, ny, nz);
mesh = enrichMesh3D(mesh,2); % Mesh detail level = 1 (for IFE).
mesh = genIntfMesh3D(mesh,pde.intf);
fem = genNedFEM3D(mesh,bc);

%% 1. use vectorization to generate face and face2elem data structure

% reoder the nodes such that the nodes of interface elements appear at the
% end

IntMeshNode = unique(reshape(mesh.t(mesh.tLoc<0,:),[],1));
MeshNodeID = true(size(mesh.p,1),1);
MeshNodeID(IntMeshNode) = false;
NodeOrderNew = zeros(size(mesh.p,1),1);
NodeOrderNew(MeshNodeID) = 1:sum(MeshNodeID);
NodeOrderNew(~MeshNodeID) = sum(MeshNodeID)+1:size(mesh.p,1);

mesh.t = NodeOrderNew(mesh.t);
mesh.e = NodeOrderNew(mesh.e);
[~,NodeOrderNewSortID] = sort(NodeOrderNew);
mesh.p = mesh.p(NodeOrderNewSortID,:);
mesh.pLoc = mesh.pLoc(NodeOrderNewSortID);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate edge structure of elements
% reorder the DoFs: put the DoFs of the edges of non-interface elements (sans those of interface elements)
% at the begining, then cancate the new edges including interface edges and
% non-interface edges

IntNodeID = [unique(reshape(mesh.t(mesh.tLoc<0,:),[],1));...
    (size(mesh.p,1)+1:size(mesh.p,1) + size(mesh.eIntP,1))'];
AllNodeID = zeros(size(node,1),1);
AllNodeID(IntNodeID) = 1:length(IntNodeID);
tfaceNew = AllNodeID(tface);
[face2EDoF,Newedge,face2EDoFSign] = dof3edge(tfaceNew);
% IMPORTANT: the Newedge contains those inside elements (4 intersection points)
% the edge orientation of dof3edge is from 1->2, 2->3, 3->1 of nodes of tfacenew
% face2EDoFSign assigns a sign such that orientation is from small to large 
Newedge = IntNodeID(Newedge); % transfer the node index of edges to the global ones
% the new index of non-interface edges of interface elements (starting from 1)
NewedgeNintID = find(abs(sum(vSign(Newedge),2))==2);

OldedgeID = unique(reshape(mesh.t_e(mesh.tLoc<0,:),[],1));
Oldedge = mesh.e(OldedgeID,:);
OldEtmp = abs(sum(vSign(Oldedge),2))==2;
OldedgeNintID = OldedgeID(OldEtmp); 
% the old index of non-interface edges of interface elements
OldedgeNint = Oldedge(OldEtmp,:);
[OldedgeNintSort,OldedgeNintSortID] = sortrows(OldedgeNint);
% ~ should be exactly the same as Newedge(NewedgeNintID)
% therefore OldedgeNintSortID should be the new index of all the
% non-interface edges of interface elements
BasicEdge = ones(size(mesh.e,1),1); BasicEdge(OldedgeID) = false;
BasicEdge = logical(BasicEdge);
% basically the non-interface edges sans edges of interface elements
BasicEdgeNum = sum(BasicEdge);
TotalOldEdge = zeros(size(mesh.e,1),1);
TotalOldEdge(BasicEdge) = 1:BasicEdgeNum;
TotalOldEdge(OldedgeNintID(OldedgeNintSortID)) = BasicEdgeNum + NewedgeNintID; 
% map OldedgeNintSortID to global ones and put them at the location of
% the original non-interface edges
% NOTE: at the location of interface edges, TotalOldEdge still contains zero

% generate new global edge DoF structure
NumEdge = size(mesh.e,1) - size(Oldedge,1) + size(Newedge,1);
% NumEdge should equal size(mesh.e,1) + sum(mesh.eLoc<0) + 2*sum(mesh.fLoc<0)
gdof = zeros(NumEdge,2); 
gdof(1:BasicEdgeNum,:) = mesh.e(BasicEdge,:);
gdof(BasicEdgeNum+1:end,:) = Newedge;
NintElem = mesh.tLoc>0;
g2ldofNint = TotalOldEdge(mesh.t_e(NintElem,:));
face2EDoF = face2EDoF + BasicEdgeNum;
% each row contains the DoFs of each element

% generate orientation structure for each edge
t_e_orit = zeros(size(g2ldofNint));
e_ind = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
for i = 1:6
    d = mesh.p(mesh.t(NintElem,e_ind(i,2)),:) - mesh.p(mesh.t(NintElem,e_ind(i,1)),:);
    %d_crect = mesh.p(mesh.e(mesh.t_e(NintElem,i),2),:) - mesh.p(mesh.e(mesh.t_e(NintElem,i),1),:);
    d_crect = node(gdof(g2ldofNint(:,i),2),:) - node(gdof(g2ldofNint(:,i),1),:);
    t_e_orit(:,i) = sign(sum(d.*d_crect,2));
end

meshI.tface = tface;
meshI.tface2elem = tface2elem;
meshI.iface = iface;
meshI.iface2elem = iface2elem;
meshI.face2elemLoc = face2elemLoc;
meshI.PolyVolume = PolyVolume;
% meshI.vSign = vSign;
meshI.node = node;
meshI.tetElem = tetElem;
meshI.tetElemLoc = tetElemLoc;
meshI.idx2cube = idx2cube;
meshI.tetVolume = volume;
meshI.face2EDoF = face2EDoF;
meshI.face2EDoFSign = face2EDoFSign;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.use face, face2elem and face2EDoF data structure to compute ife functions

face2elem = meshI.tface2elem;
tface = meshI.tface;
N = size(node, 1); NEdof = NumEdge;
NF = size(tface,1); NC = size(mesh.t, 1);

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
poly2node = sparse(face2elem(:)*ones(1,3), tface(:), 1, NP,N);
poly2node = (poly2node > 0);
NV = poly2node*ones(N, 1);

poly2edge = sparse(face2elem(:)*ones(1,3), face2EDoF(:), 1, NP,NEdof);
poly2edge = (poly2edge > 0);
NE = poly2edge*ones(NEdof, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate local IFE basis functions which are not associated with any DoFs
% IFEbasis(:,i,:,:) for minus subdomain (i=1) or plus domain (i=2)
iface = meshI.iface;
iface2elem = meshI.iface2elem;
[iface2elemReduce, uniID] = unique(iface2elem);
% length(uniID) should be NP
ifaceReduce = iface(uniID,:);
[IntFnormal,IntFarea,unitNormal,unittgt1,unittgt2] = facenormaltgt(node,ifaceReduce);
[Fnormal,Farea,FunitNormal,Funittgt1,Funittgt2] = facenormaltgt(node,tface);
% IFEbasis contains the associated IFE function on each face
% IFEbasisA is for curl term and IFEbasisB is for u term
IFEbasisA = zeros(size(face2elem,1),3,3);
% for two pieces of an ife function, only the normal vector is different
% on the minus subdomain (1) ife = n; on the plus subdomain (2) ife = n*bm/bp
IFEbasisA(:,1,:) = unitNormal(face2elem,:);
IFEbasisA(:,2,:) = unittgt1(face2elem,:);
IFEbasisA(:,3,:) = unittgt2(face2elem,:);
IFEbasisB = IFEbasisA;
piece1tmp = (face2elemLoc==1);
IFEbasisA(~piece1tmp,2,:) = IFEbasisA(~piece1tmp,2,:)*am/ap;
IFEbasisA(~piece1tmp,3,:) = IFEbasisA(~piece1tmp,3,:)*am/ap;
IFEbasisB(~piece1tmp,1,:) = IFEbasisB(~piece1tmp,1,:)*bm/bp;
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
BIFEB =  zeros(length(face2elem),3,3);
% BIFEA contains  (alpha*(grad vi).(xf -xm))for i=1,2,3 (vi is the test function)
BIFEA(:,1) = sum(squeeze(IFEbasisA(:,1,:)).*(Xf - Xmf),2);
BIFEA(:,2) = sum(squeeze(IFEbasisA(:,2,:)).*(Xf - Xmf),2);
BIFEA(:,3) = sum(squeeze(IFEbasisA(:,3,:)).*(Xf - Xmf),2);
% no need multiply area for this case as computing curl generates |F|^(-1)
BIFEA(piece1tmp,:) = am*BIFEA(piece1tmp,:);
BIFEA(~piece1tmp,:) = ap*BIFEA(~piece1tmp,:);
% BIFEB contains  wh=(beta*(curl vi)times(xf -xm))times n for i=1,2,3
% position to match the structure of mytimes
BIFEB(:,:,1) = -cross(cross(squeeze(IFEbasisB(:,1,:)),(Xf - Xmf)),FunitNormal)/2;
BIFEB(:,:,2) = -cross(cross(squeeze(IFEbasisB(:,2,:)),(Xf - Xmf)),FunitNormal)/2;
BIFEB(:,:,3) = -cross(cross(squeeze(IFEbasisB(:,3,:)),(Xf - Xmf)),FunitNormal)/2;
% devided by 2 is because curl(c times x) = 2c
BIFEB(piece1tmp,:) = bm*BIFEB(piece1tmp,:);
BIFEB(~piece1tmp,:) = bp*BIFEB(~piece1tmp,:);
% project the 3D vector onto each face
FNT = zeros(length(face2elem),3,3);
FNT(:,1,:) = Funittgt1;
FNT(:,2,:) = Funittgt2;
FNT(:,3,:) = FunitNormal;
BIFEB = mytimes(FNT,BIFEB); 
BIFEB = BIFEB(:,1:2,:);
%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%
% for BIFEB (u):
VdiagTmp1 = (bm*PolyVolume(:,1) + bm^2/bp*PolyVolume(:,2)).^(-1);
VdiagTmp2 = (bm*PolyVolume(:,1) + bp*PolyVolume(:,2)).^(-1);
VdiagTmp3 = (bm*PolyVolume(:,1) + bp*PolyVolume(:,2)).^(-1);
Vdiag1 = VdiagTmp1(face2elem);
Vdiag2 = VdiagTmp2(face2elem);
Vdiag3 = VdiagTmp3(face2elem);
Mdiag = zeros(size(Vdiag1,1),3,3);
Mdiag(:,1,1) = Vdiag1; Mdiag(:,2,2) = Vdiag2; Mdiag(:,3,3) = Vdiag3;
%%%%%%%%%%%%%%%%
% compute [-(ym-y1),(xm-x1)]
VecXm2Pt = (Xf - node(tface(:,1),:))/2;
VecXm2Pt = mytimes(FNT,VecXm2Pt);
VecXm2Pt = VecXm2Pt(:,[2,1]); VecXm2Pt(:,1) = -VecXm2Pt(:,1);
% compute [-(ym-y1),(xm-x1)].wh
BIFEB1 = zeros(size(VecXm2Pt));
for i = 1:3
    BIFEB1(:,i) = sum(VecXm2Pt.*squeeze(BIFEB(:,:,i)),2);
end
BIFEB1 = mytimes(Mdiag,BIFEB1);
% compute inv([t1;t2]) where t1 and t2 are for the 1st and 3rd edges
MatT = zeros(length(face2elem),3,2);
MatT(:,:,1) = node(tface(:,2),:) - node(tface(:,1),:);
MatT(:,:,2) = node(tface(:,1),:) - node(tface(:,3),:);
MatT = mytimes(FNT,MatT); 
MatT = MatT(:,1:2,:); MatTInv = zeros(size(MatT));
MatTInv(:,1,1) = MatT(:,2,2)./(-2*Farea); % generate inverse matrix
MatTInv(:,2,2) = MatT(:,1,1)./(-2*Farea); 
MatTInv(:,1,2) = -MatT(:,1,2)./(-2*Farea); 
MatTInv(:,2,1) = -MatT(:,2,1)./(-2*Farea); 
clear MatT
BIFEB2 = mytimes(MatTInv,BIFEB,[2,2]);  % T^(-1)w_h 
BIFEB2 = mytimes(Mdiag,permute(BIFEB2,[1,3,2]));
BIFEB2 = permute(BIFEB2,[1,3,2]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute projections to IFE spaces and form the matrices
face2EDoF = meshI.face2EDoF;
face2EDoFSign = meshI.face2EDoFSign;
IFNT = zeros(length(ifaceReduce),3,3);
IFNT(:,:,1) = unitNormal;
IFNT(:,:,2) = unittgt1;
IFNT(:,:,3) = unittgt2;

nnz = sum(NE.^2);
iiP = zeros(nnz, 1);
jjP = zeros(nnz, 1);
ssKC = zeros(nnz, 1);
ssMU = zeros(nnz, 1);
indexP = 0;
iiF = zeros(nnz, 1);
jjF = zeros(nnz, 1);
ssSC = zeros(nnz, 1);
ssSU = zeros(nnz, 1);
indexF = 0;
une = unique(NE);
hF = h(face2elem);
b = zeros(NEdof, 1);

for kk = 1:length(une) % group polys according to their # of nodes
    
    tic
    %% generate IFE functions including the IFE basis functions on each interface
    %% elements which are not associated with any DoFs, and projections of gradients
    %% and the IFE functions with matching face average
    ne = une(kk);
    % some data related to element
    isCurrentFace = (NE(face2elem) == ne);
    CurrentFace = tface(isCurrentFace,:);
    CurrentFaceE = face2EDoF(isCurrentFace,:);
    CurrentFaceS = face2EDoFSign(isCurrentFace,:);
    CNF = size(CurrentFaceE,1);
    Currentfaceloc = face2elemLoc(isCurrentFace);
    Currentface2elem = face2elem(isCurrentFace);
    CurrentFnormal = FunitNormal(isCurrentFace,:);
    FareaCurrent = Farea(isCurrentFace);
    % some data related to face
    isCurrentElem = false(NP, 1);
    isCurrentElem(face2elem(isCurrentFace)) = true;
    CurrentPolyVolume = PolyVolume(isCurrentElem,:);
    IFNTCurrent = IFNT(isCurrentElem,:,:);
    % Current elem to node matrix
    currentPoly2edge = poly2edge(isCurrentElem, :);
    CNP = size(currentPoly2edge, 1);
    currentPolyLocalIdx = zeros(NP, 1);
    currentPolyLocalIdx(isCurrentElem) = 1:CNP;
    
    % currentElem(i,:) contains the edge index of the (current) i-th element
    [I, ~] = find(currentPoly2edge');
    currentElem = reshape(I, ne, [])';
    % localIdx(i,j) contains the (current) i-th element having the
    % node(dof=1,2,3,4) on the j-th location (j is the global node index)
    localIdx = sparse(repmat((1:CNP)', 1, ne), currentElem, ones(CNP, 1)*(1:ne), CNP, NEdof);
    clear currentPoly2edge
    
    %% Deal with current triangle face cases
    %% (1) for projection of curl u
    subs1 = currentPolyLocalIdx(face2elem(isCurrentFace));
    subs2 = [ones(CNF, 1); 2*ones(CNF, 1); 3*ones(CNF, 1)];
    val = BIFEA(isCurrentFace,:);
    
    % compute the projection of gradients
    GPA1 = zeros(CNP,3,ne); 
    % coefficient of n t1 and t2 for the subelement on the subdomain 1 
    for m = 1:3
        subs3 = full(localIdx((CurrentFaceE(:, m) - 1)*CNP + subs1));
        SignCurrentE = repmat(CurrentFaceS(:,m),3,1);
        % subs1: new order (range in 1:CNP) of the element index
        % subs2: three components in the vector
        % subs3: the local index (DoF) (w.r.t. element) of the m-th edge of this face
        GPA1 = GPA1 + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)],...
            SignCurrentE.*val(:), [CNP, 3, ne]);
    end
    GPA2 = GPA1; GPA2(:,2,:) = am/ap*GPA2(:,2,:); GPA2(:,3,:) = am/ap*GPA2(:,3,:);
    GPA1 = mytimes(IFNTCurrent,GPA1);
    GPA2 = mytimes(IFNTCurrent,GPA2);
    % the updated GP contains three component values of the projection vector   
    % generate IFE vectors on each face
    % clear currentPolyLocalIdx localIdx
    % clear subs1 subs2 subs3
    
    % generate the data for computing stabilization matrices
    CurrentFace2DoF = zeros(NP,ne);
    CurrentFace2DoF(isCurrentElem,:) = currentElem;
    CurrentFace2DoF = CurrentFace2DoF(Currentface2elem,:);
    % CurrentFace2DoF contains the DoFs of the element associated with each face
    
    % compute curl uh.n on each face (uh is an IFE function)
    GPA1tmp = zeros(NP,3,ne); GPA2tmp = zeros(NP,3,ne);
    GPA1tmp(isCurrentElem,:,:) = GPA1; GPA2tmp(isCurrentElem,:,:) = GPA2;
    GAface = zeros(CNF,3,ne);
    GAface(Currentfaceloc==1,:,:) = GPA1tmp(Currentface2elem(Currentfaceloc==1),:,:);
    GAface(Currentfaceloc==2,:,:) = GPA2tmp(Currentface2elem(Currentfaceloc==2),:,:);
    clear GPA1tmp GPA2tmp     
    GAfaceP = zeros(size(GAface,1),ne);
    for n = 1:ne
        GAfaceP(:,n) =  sum(squeeze(GAface(:,:,n)).*CurrentFnormal,2);
    end
    
    % compute rot_F uh on each face according to DoFs
    TrueCurlTmp = CurrentFaceS./FareaCurrent;
    TrueCurl = zeros(CNF,ne);
    for i = 1:ne
        for j = 1:3
            IDij = (CurrentFaceE(:,j) == CurrentFace2DoF(:,i));
            TrueCurl(IDij,i) = TrueCurlTmp(IDij,j);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (2) for projection of u
    val = BIFEB1(isCurrentFace,:);   
    GPB1 = zeros(CNP,3,ne);
    % coefficient of n t1 and t2 for the subelement on the subdomain 1
    for m = 1:3
        subs3 = full(localIdx((CurrentFaceE(:, m) - 1)*CNP + subs1));
        SignCurrentE = repmat(CurrentFaceS(:,m),3,1);
        % subs1: new order (range in 1:CNP) of the element index
        % subs2: three components in the vector
        % subs3: the local index (DoF) (w.r.t. element) of the m-th edge of this face
        GPB1 = GPB1 + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)],...
            SignCurrentE.*val(:), [CNP, 3, ne]);
    end
    
    edgetmp = [1,3];
    for mm = 1:2
        val = squeeze(BIFEB2(isCurrentFace,mm,:));
        m = edgetmp(mm);
        subs3 = full(localIdx((CurrentFaceE(:, m) - 1)*CNP + subs1));
        SignCurrentE = repmat(CurrentFaceS(:,m),3,1);
        Areatmp = repmat(FareaCurrent,3,1);
        % subs1: new order (range in 1:CNP) of the element index
        % subs2: three components in the vector
        % subs3: the local index (DoF) (w.r.t. element) of the m-th edge of this face
        GPB1 = GPB1 + accumarray([repmat(subs1,3,1), subs2, repmat(subs3, 3, 1)],...
            SignCurrentE.*val(:).*Areatmp, [CNP, 3, ne]);
    end    
    GPB2 = GPB1; GPB2(:,2,:) = am/ap*GPB2(:,2,:); GPB2(:,3,:) = am/ap*GPB2(:,3,:);
    GPB1 = mytimes(IFNTCurrent,GPB1);
    GPB2 = mytimes(IFNTCurrent,GPB2);
    % the updated GP contains three component values of the projection vector   
    % generate IFE vectors on each face
    clear currentPolyLocalIdx localIdx
    clear subs1 subs2 subs3
    
    % compute uh^(\tau) (projection onto each face)
    GPB1tmp = zeros(NP,3,ne); GPB2tmp = zeros(NP,3,ne);
    GPB1tmp(isCurrentElem,:,:) = GPB1; GPB2tmp(isCurrentElem,:,:) = GPB2;
    GBface = zeros(CNF,3,ne);
    GBface(Currentfaceloc==1,:,:) = GPB1tmp(Currentface2elem(Currentfaceloc==1),:,:);
    GBface(Currentfaceloc==2,:,:) = GPB2tmp(Currentface2elem(Currentfaceloc==2),:,:);
    clear GPA1tmp GPA2tmp     
    GBfaceP = zeros(size(GBface));
    for n = 1:ne
        GBfaceP(:,:,n) = GBface(:,:,n) - sum(GBface(:,:,n).*CurrentFnormal,2).*CurrentFnormal;
    end
    clear GBface
    GBfaceP = mytimes(FNT(isCurrentFace,:,:),GBfaceP);
    GBfaceP = GBfaceP(:,1:2,:);
    
    % compute u^(tau)(xm) uh on each face according to DoFs
    TrueUtTmp = zeros(CNF,2,ne);
    TrueUt = zeros(CNF,2,ne);
    MatTInvCurrent = MatTInv(isCurrentFace,:,:);
    for i = 1:ne
        TrueUt(:,:,i) = TrueUt(:,:,i) + TrueCurl(:,i).*VecXm2Pt(isCurrentFace,:);
        for j = 1:2
            jj = edgetmp(j);
            IDij = (CurrentFaceE(:,jj) == CurrentFace2DoF(:,i));
            TrueUtTmp(IDij,:,i) = MatTInvCurrent(IDij,j,:).*CurrentFaceS(IDij,jj);
        end
    end
    TrueUt = TrueUt + TrueUtTmp;
    clear TrueUtTmp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate stiffness and mass matrices
    for n = 1:ne
        for m = 1:ne
            iiP(indexP+1:indexP + CNP) = currentElem(:, n);
            jjP(indexP+1:indexP + CNP) = currentElem(:, m);
            ssKC(indexP+1:indexP + CNP) = am*dot(GPA1(:,:,n), GPA1(:,:,m),2).*CurrentPolyVolume(:,1)+...
                ap*dot(GPA2(:,:,n), GPA2(:,:,m),2).*CurrentPolyVolume(:,2);
            ssMU(indexP+1:indexP + CNP) = bm*dot(GPB1(:,:,n), GPB1(:,:,m),2).*CurrentPolyVolume(:,1)+...
                bp*dot(GPB2(:,:,n), GPB2(:,:,m),2).*CurrentPolyVolume(:,2);
            indexP = indexP + CNP;
            
            iiF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, n);
            jjF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, m);
            ssSC(indexF+1:indexF + CNF) = dot((TrueCurl(:,m) - GAfaceP(:,m)),...
                (TrueCurl(:,n) - GAfaceP(:,n)),2).*...
                Farea(isCurrentFace).*hF(isCurrentFace);
            ssSU(indexF+1:indexF + CNF) = ( dot((TrueUt(:,1,m) - GBfaceP(:,1,m)),...
                (TrueUt(:,1,n) - GBfaceP(:,1,n)),2) + dot((TrueUt(:,2,m) - GBfaceP(:,2,m)),...
                (TrueUt(:,2,n) - GBfaceP(:,2,n)),2) ).*Farea(isCurrentFace).*hF(isCurrentFace);
            indexF = indexF + CNF;
        end
    end
    
    %% form the right-hand side vector
    idx2cubeNew = -mesh.tLoc(idx2cube);
    TetID = isCurrentElem(idx2cubeNew);
    currentvolume = volume(TetID);
    %%% DoFs of elements for tet
    currentTetDoF = zeros(size(isCurrentElem,1),ne);
    currentTetDoF(isCurrentElem,:) = currentElem;
    currentTetDoF = currentTetDoF(idx2cubeNew(TetID),:);
    %%% quadrature points
    X1 = node(tetElem(TetID,1),:); 
    X2 = node(tetElem(TetID,2),:); 
    X3 = node(tetElem(TetID,3),:);
    X4 = node(tetElem(TetID,4),:);
    ng = 1;
    [gx, gy, gz] = gaussPtetra(X1, X2, X3, X4, ng);
    gw = gaussWtetra(ng);
    Xm = (node(ifaceReduce(isCurrentElem,1),:) + node(ifaceReduce(isCurrentElem,2),:) +...
        node(ifaceReduce(isCurrentElem,3),:))/3;
    Xmtet = zeros(size(isCurrentElem,1),3);
    Xmtet(isCurrentElem,:) = Xm;
    Xmtet = Xmtet(idx2cubeNew(TetID),:);
    %%% projection of u
    Bas1 = zeros(size(isCurrentElem,1),3,ne); Bas2 = Bas1;
    Bas1(isCurrentElem,:,:) = GPB1; 
    Bas2(isCurrentElem,:,:) = GPB2; 
    Bas1 = Bas1(idx2cubeNew(TetID),:,:); Bas2 = Bas2(idx2cubeNew(TetID),:,:);
    Bas = zeros(sum(TetID),3,ne);
    piecetmp = tetElemLoc(TetID);
    Bas(piecetmp==1,:,:) = Bas1(piecetmp==1,:,:);
    Bas(piecetmp==2,:,:) = Bas2(piecetmp==2,:,:);
    
    
    ft1 = pde.f1(gx,gy,gz);
    ft2 = pde.f2(gx,gy,gz);
    ft3 = pde.f3(gx,gy,gz);
    Currentb = zeros(size(ft1));
    cnt = size(ft1,1);
    bii = ne*cnt;
    count = 0 ;
    for i = 1:ne
        Basi = squeeze(Bas(:,:,i));
        Currentb(:,i) = sum(gw*(ft1.*Basi(:,1)+ft2.*Basi(:,2)+ft3.*Basi(:,3)),2).*...
            currentvolume;
        bii(count+1:count+cnt) = currentTetDoF(:,i);
        count = count+cnt;
    end
    b = b + sparse(bii,1,reshape(Currentb,[],1),NEdof,1);
    
    toc    
    
end

Kcurl = sparse(iiP, jjP, ssKC, NEdof, NEdof);
Scurl = sparse(iiF, jjF, ssSC, NEdof, NEdof);
MU = sparse(iiP, jjP, ssMU, NEdof, NEdof);
SU = sparse(iiF, jjF, ssSU, NEdof, NEdof);

MatI = Kcurl + Scurl + MU +SU;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate matrix on non-interface elements (for test)
feEvalBas1 = @EvalNed1Bas3D;
feEvalBas2 = @EvalNed1Bas3D;

dof1 = 6; dof2 = 6; nloc = dof1*dof2; 
ntID = find(mesh.tLoc > 0); ntN = length(ntID);

%% 1. Matrix on noninterface elements
MatInd = 1;
AN = fem.area(ntID); 
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); gw = fem.gw;
XN = zeros(nloc*ntN, 1);

coefN = feval(pde.A,gxN,gyN,gzN);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);

for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 1).*t_e_orit(:,i);
    Ibasy{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 2).*t_e_orit(:,i);
    Ibasz{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 3).*t_e_orit(:,i);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 1).*t_e_orit(:,j);
    Jbasy{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 2).*t_e_orit(:,j);
    Jbasz{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 3).*t_e_orit(:,j);
end

IN = reshape(repmat(g2ldofNint(:,1:6),6,1),nloc*ntN,1);
JN = repmat(reshape(g2ldofNint(:,1:6),dof2*ntN,1),6,1);
ind = 0;
for i = 1:dof1
    for j = 1:dof2
        XN(ind+1:ind+ntN) = AN.*(sum(((Ibasx{i}.*(coefN.*Jbasx{j})).*gw'),2) + ...
            sum(((Ibasy{i}.*(coefN.*Jbasy{j})).*gw'),2) + ...
            sum(((Ibasz{i}.*(coefN.*Jbasz{j})).*gw'),2));
        ind = ind + ntN;
    end
end
ID = find(XN~=0); 
SN = sparse(IN(ID),JN(ID),XN(ID),NEdof,NEdof);

StiffCurl = SN + Kcurl + Scurl;

%%%%%%

MatInd = 0;
AN = fem.area(ntID); 
gxN = fem.gx(ntID,:); gyN = fem.gy(ntID,:); gzN = fem.gz(ntID,:); gw = fem.gw;
XN = zeros(nloc*ntN, 1);

coefN = feval(pde.A,gxN,gyN,gzN);
Ibasx = cell(dof1,1); Ibasy = cell(dof1,1); Ibasz = cell(dof1,1); 
Jbasx = cell(dof2,1); Jbasy = cell(dof2,1); Jbasz = cell(dof2,1);

for i = 1:dof1
    Ibasx{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 1).*t_e_orit(:,i);
    Ibasy{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 2).*t_e_orit(:,i);
    Ibasz{i} = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, MatInd, 3).*t_e_orit(:,i);
end
for j = 1:dof2
    Jbasx{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 1).*t_e_orit(:,j);
    Jbasy{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 2).*t_e_orit(:,j);
    Jbasz{j} = feEvalBas2(fem.bas, ntID, gxN, gyN, gzN, j, MatInd, 3).*t_e_orit(:,j);
end

IN = reshape(repmat(g2ldofNint(:,1:6),6,1),nloc*ntN,1);
JN = repmat(reshape(g2ldofNint(:,1:6),dof2*ntN,1),6,1);
ind = 0;
for i = 1:dof1
    for j = 1:dof2
        XN(ind+1:ind+ntN) = AN.*(sum(((Ibasx{i}.*(coefN.*Jbasx{j})).*gw'),2) + ...
            sum(((Ibasy{i}.*(coefN.*Jbasy{j})).*gw'),2) + ...
            sum(((Ibasz{i}.*(coefN.*Jbasz{j})).*gw'),2));
        ind = ind + ntN;
    end
end
ID = find(XN~=0); 
MN = sparse(IN(ID),JN(ID),XN(ID),NEdof,NEdof);

MassU = MN + MU + SU;

%%%%%%%%%%

dof = 6;  nloc = dof;
X = zeros(nloc*ntN, 1);

fN1 = feval(pde.exactu1,gxN,gyN,gzN);
fN2 = feval(pde.exactu2,gxN,gyN,gzN);
fN3 = feval(pde.exactu3,gxN,gyN,gzN);
ind = 0;
I = reshape(g2ldofNint(:,1:6),nloc*ntN,1);
for i = 1:dof
    ibas1 = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, 0, 1);
    ibas2 = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, 0, 2);
    ibas3 = feEvalBas1(fem.bas, ntID, gxN, gyN, gzN, i, 0, 3);
    X(ind+1:ind+ntN) = AN.*sum((ibas1.*fN1+ibas2.*fN2+ibas3.*fN3).*gw',2).*fem.t_e_orit(ntID,i);
    ind = ind + ntN;
end
rhsN = sparse(I,1,X,NEdof,1);

rhs = rhsN + b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gew,gex,gey,gez] = gaussPedge(node(gdof(:,1),:),node(gdof(:,2),:),1);


tu1 = sum(feval(pde.exactu1,gex,gey,gez).*gew,2);
tu2 = sum(feval(pde.exactu2,gex,gey,gez).*gew,2);
tu3 = sum(feval(pde.exactu3,gex,gey,gez).*gew,2);

tgt = node(gdof(:,2),:) - node(gdof(:,1),:);
tgt = tgt./sum(tgt.^2,2).^(1/2);
tu = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);

Atotal = StiffCurl + MassU;
%b = ones(size(MassU,1),1);
x0 = zeros(size(b));

bcind = fem.bcind;
[bc,mapper] = boundaryEdge3D(node,gdof,bcind);
bdidx = zeros(NEdof,1); 
isBdEdge = true(NEdof,1);
isBdEdge(mapper) = false;
bdidx(isBdEdge) = 1;
Tbd = spdiags(bdidx,0,NEdof,NEdof);
T = spdiags(1-bdidx,0,NEdof,NEdof);
A = T*Atotal*T + Tbd;

ub = tu;
ub(mapper) = 0;
rhsB = Atotal*ub;
f = rhs - rhsB;
f(isBdEdge) = tu(isBdEdge);

option.outsolver = 'cg';
alpha = am*ones(size(gdof,1),1);
% alpha(femI.eLoc>0) = ap;
% alpha(femI.eLoc==0) = (am+ap)/2;
beta = bm*ones(size(gdof,1),1);
% beta(femI.eLoc>0) = bp;
% beta(femI.eLoc==0) = (bm+bp)/2;
option.alpha = alpha;
option.beta = beta;
option.solver = 'amg';
edge = gdof;
option.isBdEdge = isBdEdge;
option.smoother = 'BD';
D = diag(A); D = D(1:BasicEdgeNum);
M = sparse(1:BasicEdgeNum,1:BasicEdgeNum,D,NEdof,NEdof);
M(BasicEdgeNum+1:end,BasicEdgeNum+1:end) = M(BasicEdgeNum+1:end,BasicEdgeNum+1:end) +...
    A(BasicEdgeNum+1:end,BasicEdgeNum+1:end);
option.M = M;
%option.M = A(BasicEdgeNum+1:end,BasicEdgeNum+1:end);
option.blkId = BasicEdgeNum;
[x,info] = amgMaxwell(A,f,node,edge,option);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate exact solution 

% [gew,gex,gey,gez] = gaussPedge(node(gdof(:,1),:),node(gdof(:,2),:),1);
% 
% 
% tu1 = sum(feval(pde.exactu1,gex,gey,gez).*gew,2);
% tu2 = sum(feval(pde.exactu2,gex,gey,gez).*gew,2);
% tu3 = sum(feval(pde.exactu3,gex,gey,gez).*gew,2);
% 
% tgt = node(gdof(:,2),:) - node(gdof(:,1),:);
% tgt = tgt./sum(tgt.^2,2).^(1/2);
% tu = tu1.*tgt(:,1) + tu2.*tgt(:,2) + tu3.*tgt(:,3);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% test
% DC = full(StiffCurl*tu);
% DU = full(SU*tu);


% for i = 1:CNP
%     vv = sum(tu(currentElem(1,:)).*(squeeze(GPA1(1,:,:))'),1);
%     if sum(abs(vv)) > 10^(-8)
%         stp
%     end
% end