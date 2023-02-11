function femI = genNedVIFEM3D(meshI,mesh,pde,option)
%% use face and face2elem data structure to compute ife functions
bm = pde.bm; bp = pde.bp;
%bm = 1; bp = 1;
am = pde.am; ap = pde.ap;
node = [mesh.p;mesh.eIntP];
face2elem = meshI.tface2elem;
iface = meshI.iface;
tface = meshI.tface;
iface2elem = meshI.iface2elem;
face2elemLoc = meshI.face2elemLoc;
PolyVolume = meshI.PolyVolume;
idx2cube = meshI.idx2cube;
volume = meshI.tetVolume;
tetElem = meshI.tetElem;
tetElemLoc = meshI.tetElemLoc;

h = (PolyVolume(:,1) + PolyVolume(:,2)).^(1/3);
vSign = [mesh.pLoc;zeros(size(mesh.eIntP,1),1)];

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
% t_e_orit = zeros(size(g2ldofNint));
% e_ind = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
% for i = 1:6
%     d = mesh.p(mesh.t(NintElem,e_ind(i,2)),:) - mesh.p(mesh.t(NintElem,e_ind(i,1)),:);
%     %d_crect = mesh.p(mesh.e(mesh.t_e(NintElem,i),2),:) - mesh.p(mesh.e(mesh.t_e(NintElem,i),1),:);
%     d_crect = node(gdof(g2ldofNint(:,i),2),:) - node(gdof(g2ldofNint(:,i),1),:);
%     t_e_orit(:,i) = sign(sum(d.*d_crect,2));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NEdof = NumEdge;
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
clear isCutPoly cutPolyIdx

% prepare new data structure to compute projections
poly2edge = sparse(face2elem(:)*ones(1,3), face2EDoF(:), 1, NP,NEdof);
poly2edge = (poly2edge > 0);
NE = poly2edge*ones(NEdof, 1);
%centroid = poly2node*node./[NV, NV, NV];
%KP = pde.A(centroid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate local IFE basis functions which are not associated with any DoFs
% IFEbasis(:,i,:,:) for minus subdomain (i=1) or plus domain (i=2)

[~, uniID] = unique(iface2elem);
% length(uniID) should be NP
ifaceReduce = iface(uniID,:);
[~,~,unitNormal,unittgt1,unittgt2] = facenormaltgt(node,ifaceReduce);
[~,Farea,FunitNormal,Funittgt1,Funittgt2] = facenormaltgt(node,tface);
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

ElemID = cell(1,2);
BasA1 = cell(1,2); BasA2 = cell(1,2);
BasB1 = cell(1,2); BasB2 = cell(1,2);
ElemDoF = cell(1,2);

for kk = 1:length(une) % group polys according to their # of nodes
    
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
    GPB2 = GPB1;
    GPB2(:,1,:) = bm/bp*GPB1(:,1,:); %GPB2(:,2,:) = GPB2(:,2,:); GPB2(:,3,:) = GPB2(:,3,:);
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
            
            if isfield(option,'w1')
                w1 = option.w1;
                w2 = option.w2;
            else
                w1 = 1; w2 =1;
            end
            iiF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, n);
            jjF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, m);
            ssSC(indexF+1:indexF + CNF) = dot((TrueCurl(:,m) - GAfaceP(:,m)),...
                (TrueCurl(:,n) - GAfaceP(:,n)),2).*...
                Farea(isCurrentFace)*w1.*hF(isCurrentFace); % 20 is better for large jump
            ssSU(indexF+1:indexF + CNF) = ( dot((TrueUt(:,1,m) - GBfaceP(:,1,m)),...
                (TrueUt(:,1,n) - GBfaceP(:,1,n)),2) + dot((TrueUt(:,2,m) - GBfaceP(:,2,m)),...
                (TrueUt(:,2,n) - GBfaceP(:,2,n)),2) ).*Farea(isCurrentFace)*w2;%.*hF(isCurrentFace);
                % weight 0.1 is good for time harmonic equations
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
    
    ElemID{kk} = isCurrentElem;
    BasA1{kk} = GPA1; BasA2{kk} = GPA2;
    BasB1{kk} = GPB1; BasB2{kk} = GPB2;
    ElemDoF{kk} = currentElem;
    
end

Kcurl = sparse(iiP, jjP, ssKC, NEdof, NEdof);
Scurl = sparse(iiF, jjF, ssSC, NEdof, NEdof);
MU = sparse(iiP, jjP, ssMU, NEdof, NEdof);
SU = sparse(iiF, jjF, ssSU, NEdof, NEdof);

KI = Kcurl + Scurl; 
SI = MU + SU;

eLoc = zeros(size(gdof,1),1);
eLoc(min(vSign(gdof),[],2)==1) = 1;
eLoc(max(vSign(gdof),[],2)==-1) = -1;


femI.BasA1 = BasA1; femI.BasA2 = BasA2;
femI.BasB1 = BasB1; femI.BasB2 = BasB2;
femI.KI = KI; 
femI.SI = SI; 
femI.b = b;
femI.ElemID = ElemID;
femI.ElemDoF = ElemDoF;
%femI.t_e_orit = t_e_orit;
femI.g2ldofNint = g2ldofNint;
femI.gdof = gdof;
femI.BasicEdgeNum = BasicEdgeNum;
femI.eLoc = eLoc;
femI.une = une;