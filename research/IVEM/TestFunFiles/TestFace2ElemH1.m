
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. use vectorization to generate face and face2elem data structure
%% this vectorization code does not work for uniform tets

% tic
% 
% matlabTri = 1;
% N = size(mesh.p,1);
% node = [mesh.p;mesh.eIntP];
% allCutElem = mesh.t(mesh.tLoc<0,:);
% isCutElem = (mesh.tLoc<0);
% vSign = [mesh.pLoc;zeros(size(mesh.eIntP,1),1)];
% 
% isInterfaceNode = false(size(node,1), 1);
% isInterfaceNode(allCutElem(:)) = true; % include the vertices of allCutElem 
% isInterfaceNode(N+1:end) = true;    % and the cut points and the aux points
% interfaceNode = node(isInterfaceNode,:);
% matlabversion = version();
% if matlabTri == 1
%     DT = delaunayTriangulation(interfaceNode);
%     tetElem = DT.ConnectivityList;
% elseif matlabTri == 0
%     constrant1 = [mesh.e(mesh.eLoc<0,1),N-mesh.eLoc(mesh.eLoc<0)];
%     constrant2 = [N-mesh.eLoc(mesh.eLoc<0),mesh.e(mesh.eLoc<0,2)];
%     constrant = [constrant1;constrant2];
%     DT = constrainedDelaunayTetGen(interfaceNode,constrant);
%     tetElem = DT.ConnectivityList;
% end
% [tetElem, ortidx, volume] = fixorder3(interfaceNode, tetElem);
% 
% localidx2globalidx = find(isInterfaceNode); % tetElem points to interfaceNode
% tetElem = localidx2globalidx(tetElem);  % map to the global index of node
% bc = (node(tetElem(:, 1), :) + node(tetElem(:, 2), :) ...
%     + node(tetElem(:, 3), :) + node(tetElem(:, 4), :))/4.0;
% idx2cube = FindElemId(bc, mesh); % Get the cube index of each tetElem
% idx2cubeBad = (idx2cube<=0);
% idx2cube = idx2cube(~idx2cubeBad);
% tetElem = tetElem(isCutElem(idx2cube), :); % Just keep the tets in cut elem
% volume = volume(isCutElem(idx2cube), :);
% idx2cube = idx2cube(isCutElem(idx2cube));  % keep cube idx in cut elements
% % volumetet0 =  accumarray(idx2cube,volume);
% % volumetet0 = volumetet0(volumetet0>0);
% %[tetElemOld, ortidxOld, volumeOld] = fixorder3(interfaceNode, mesh.t(mesh.tLoc<0,:));
% 
% % There might exist some bad tet elements whose vertices are in
% % different cubes. We just need to keep the tet in every cut elem. So
% % here we get rid of them. (but do not handle different tets)
% NT = size(tetElem,1);
% X = reshape(node(tetElem,1), NT, 4);
% Y = reshape(node(tetElem,2), NT, 4);
% Z = reshape(node(tetElem,3), NT, 4);
% isBadTet =  (max(X,[],2) - min(X,[],2)) > 2*h - eps ... % d > 2h means
%           | (max(Y,[],2) - min(Y,[],2)) > 2*h - eps ... % in differnt cube
%           | (max(Z,[],2) - min(Z,[],2)) > 2*h - eps;
% tetElem = tetElem(~isBadTet,:);
% 
% toc

% plot to check
%showmesh3(node,tetElem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this approach is based on forlopp which is slow
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
tetElem = localidx2globalidx(tetElem);  % map to the global index of node
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try another approach for vectorization

% tic
% 
% Np = size(mesh.p,1);
% Ncube = size(mesh.t,1)/6;
% node = [mesh.p;mesh.eIntP];
% allCutElem = mesh.t(mesh.tLoc<0,:);
% isCutElem = (mesh.tLoc<0);
% vSign = [mesh.pLoc;zeros(size(mesh.eIntP,1),1)];
% intID = find(mesh.tLoc<0);
% isInterfaceNode = false(size(node,1), 1);
% isInterfaceNode(allCutElem(:)) = true; % include the vertices of allCutElem 
% isInterfaceNode(Np+1:end) = true;    % and the cut points and the aux points
% 
% for j = 1:6
%     
%     Indj = (j-1)*Ncube + 1:Ncube;
%     Indj = intersect(intID,Indj);
%     allCutElemj = mesh.t(Indj,:);
%     isInterfaceNodej = false(size(node,1), 1);
%     isInterfaceNodej(allCutElemj(:)) = true; % include the vertices of allCutElem
%     tmp = mesh.eLoc(mesh.t_e(Indj,:));
%     tmp = tmp(tmp<0);
%     %isInterfaceNodej(Np-tmp) = true;    % and the cut points and the aux points
%     interfaceNodej = node(isInterfaceNodej,:);
%     
%     DT = delaunayTriangulation(interfaceNodej);
%     tetElemj = DT.ConnectivityList;
%     [tetElemj, ortidx, volume] = fixorder3(interfaceNodej, tetElemj);
%     
%     localidx2globalidxj = find(isInterfaceNodej); % tetElem points to interfaceNode
%     tetElemj = localidx2globalidxj(tetElemj);  % map to the global index of node
%     bcj = (node(tetElemj(:, 1), :) + node(tetElemj(:, 2), :) ...
%         + node(tetElemj(:, 3), :) + node(tetElemj(:, 4), :))/4.0;
%     idx2cubej = FindElemId(bcj, mesh); % Get the cube index of each tetElem
%     idx2cubeBadj = (idx2cubej<=0);
%     idx2cubej = idx2cubej(~idx2cubeBadj);
%     tetElemj = tetElemj(~idx2cubeBadj,:);
%     [godElemj,godElemjIDa,godElemjIDb] = intersect(idx2cubej,Indj);
%     idx2cubej = idx2cubej(godElemjIDa);
%     tetElemj = tetElemj(godElemjIDa,:);
%   
% %     tetElem = tetElem(isCutElem(idx2cube), :); % Just keep the tets in cut elem
% %     volume = volume(isCutElem(idx2cube), :);
% %     idx2cube = idx2cube(isCutElem(idx2cube));  % keep cube idx in cut elements
%     
% end
% 
% toc


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

% plot to check
%trisurf(tface(c2old+1:c2,:),node(:,1),node(:,2),node(:,3))
%trisurf(iface,node(:,1),node(:,2),node(:,3))

tface((ct+1):end,:) = [];    % remove empty meomory
tface2elem((ct+1):end) = [];

face2elem = tface2elem;
face2elemLoc = [ones(c2old,1);2*ones(c2-c2old,1)]; % face2elemLoc is the same as face location
face = tface;

interfaceData.tface = tface;
interfaceData.tface2elem = tface2elem;
interfaceData.iface = iface;
interfaceData.iface2elem = iface2elem;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.use face and face2elem data structure to compute ife functions

vSign = [mesh.pLoc;zeros(size(mesh.eIntP,1),1)];
node = [mesh.p;mesh.eIntP];
face = interfaceData.tface;
face2elem = interfaceData.tface2elem;
N = size(node, 1); Ndof = N;
NF = size(face,1); NC = size(mesh.t, 1);

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


% generate local IFE basis functions which are not associated with any DoFs
% IFEbasis(:,i,:,:) for minus subdomain (i=1) or plus domain (i=2)
iface = interfaceData.iface;
tface = interfaceData.tface;
iface2elem = interfaceData.iface2elem;
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
% generate boundary integral for computing projections
BIFE =  zeros(length(face2elem),3);
% BIFE contains  (beta*grad vi.n)*|f|/3 for i=1,2,3 (vi is the test function)
BIFE(:,1) = sum(squeeze(IFEbasis(:,1,:)).*Fnormal,2)/2; 
BIFE(:,2) = sum(squeeze(IFEbasis(:,2,:)).*Fnormal,2)/2;
BIFE(:,3) = sum(squeeze(IFEbasis(:,3,:)).*Fnormal,2)/2;
% devided by 2 is because Fnormal is unitnormal*(2*area)
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute projections to IFE spaces and form the matrices
nnz = sum(NV.^2);
iiP = zeros(nnz, 1);
jjP = zeros(nnz, 1);
ssP = zeros(nnz, 1);
indexP = 0;
iiF = zeros(nnz, 1);
jjF = zeros(nnz, 1);
ssF = zeros(nnz, 1);
indexF = 0;
unv = unique(NV);
hF = h(face2elem);
b = zeros(Ndof, 1);

for kk = 1:length(unv) % group polys according to their # of nodes
    
    tic
    %% generate IFE functions including the IFE basis functions on each interface
    %% elements which are not associated with any DoFs, and projections of gradients
    %% and the IFE functions with matching face average
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
    toc
    
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
    % compute the face gradients based on DoFs
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
            
            iiF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, n);
            jjF(indexF+1:indexF + CNF) = CurrentFace2DoF(:, m);
            ssF(indexF+1:indexF + CNF) = hF(isCurrentFace).*...
                dot((gradLambda(:,:,m) - GfaceP(:,:,m)),(gradLambda(:,:,n) - GfaceP(:,:,n)),2).*...
                Farea(isCurrentFace);
            indexF = indexF + CNF;
        end
    end
    
    
    %% compute the RHS with only evaluating the integral at the centroid
    %% approach 1 (will give the optimal order but slightly increase the error)
%     currentElemSign = vSign(currentElem);
%     currentElemSign1 = currentElemSign<=0; nv1 = sum(currentElemSign1,2);
%     currentElemSign2 = currentElemSign>=0; nv2 = sum(currentElemSign2,2);
% %     nodeX = node(:,1); currentElemNodeX = nodeX(currentElem); clear nodeX
% %     nodeY = node(:,2); currentElemNodeY = nodeY(currentElem); clear nodeY
% %     nodeZ = node(:,1); currentElemNodeZ = nodeZ(currentElem); clear nodeZ
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
    
    %% approach 2
    idx2cubeNew = -mesh.tLoc(idx2cube);
    TetID = isCurrentElem(idx2cubeNew);
    currentvolume = volume(TetID);
    ElemDoF = zeros(sum(TetID),nv);
    ElemDoF(isCurrentElem,:) = currentElem; 
    ElemDoF = ElemDoF(idx2cubeNew(TetID),:);
    X1 = node(tetElem(TetID,1),:); 
    X2 = node(tetElem(TetID,2),:); 
    X3 = node(tetElem(TetID,3),:);
    X4 = node(tetElem(TetID,4),:);
    ng = 1;
    [gx, gy, gz] = gaussPtetra(X1, X2, X3, X4, ng);
    gw = gaussWtetra(ng);
    Xm = (node(ifaceReduce(isCurrentElem,1),:) + node(ifaceReduce(isCurrentElem,2),:) +...
        node(ifaceReduce(isCurrentElem,3),:))/3;
    
    Bas1 = zeros(size(isCurrentElem,1),3,nv); Bas2 = Bas1;
    Bas1(isCurrentElem,:,:) = GP1; 
    Bas2(isCurrentElem,:,:) = GP2; 
    Bas1 = Bas1(idx2cubeNew(TetID),:,:); Bas2 = Bas2(idx2cubeNew(TetID),:,:);
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
    Bas(piecetmp==1,:,:) = Bas1(piecetmp==1,:,:);
    Bas(piecetmp==2,:,:) = Bas2(piecetmp==2,:,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mass matrix
    
%     for i = 1:nv
%         for j = 1:nv
%             Basi = squeeze(Bas(:,:,i)); IFEc0i = IFEc0(:,i);
%             Basj = squeeze(Bas(:,:,j)); IFEc0j = IFEc0(:,j);
%             uhi = (gx-Xmtet(:,1)).*Basi(:,1) + (gy-Xmtet(:,2)).*Basi(:,2) +...
%                 (gz-Xmtet(:,3)).*Basi(:,3) + IFEc0i;
%             uhj = (gx-Xmtet(:,1)).*Basj(:,1) + (gy-Xmtet(:,2)).*Basj(:,2) +...
%                 (gz-Xmtet(:,3)).*Basj(:,3) + IFEc0j;
%             Currentb(:,i) = sum(gw*ft.*uhi.*uhj,2).*currentvolume;
%             bii(count+1:count+cnt) = currentTetDoF(:,i);
%             count = count+cnt;
%         end
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ft = pde.f(gx,gy,gz);
    Currentb = zeros(size(ft));
    cnt = size(ft,1);
    bii = nv*cnt;
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
    %b = b + accumarray([idx2cubeNew(TetID),currentTetDoF(:)], Currentb(:), [Ndof, 1]);
    

end
K = sparse(iiP, jjP, ssP, Ndof, Ndof);
S = sparse(iiF, jjF, ssF, Ndof, Ndof);













