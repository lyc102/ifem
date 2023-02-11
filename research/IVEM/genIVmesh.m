function meshI = genIVmesh(mesh)

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
    if i == 3
        stp = 1;
    end
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
d12 = interfaceNode(tetElem(:,2),:) - interfaceNode(tetElem(:,1),:);
d13 = interfaceNode(tetElem(:,3),:) - interfaceNode(tetElem(:,1),:);
d14 = interfaceNode(tetElem(:,4),:) - interfaceNode(tetElem(:,1),:);
d12 = sum(d12.^2,2).^(1/2); d13 = sum(d13.^2,2).^(1/2); d14 = sum(d14.^2,2).^(1/2);
ld = max([d12,d13,d14],[],2);
idPlane = (volume./ld<=10^(-16));


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


%% Get triangular faces for interrior elements and interface
localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
NT = size(tetElem,1);
tface = zeros(4*NT, 3);
tface2elem = zeros(4*NT, 1);
iface = zeros(4*NT, 3);
iface2elem = zeros(4*NT, 1);
% find the interior tet elements
isIntTet1 = min(vSign(tetElem),[], 2) == -1; % can not be == -1 as there is sliver
isIntTet2 = sum(abs(vSign(tetElem)),2) == 0;
isIntTet = isIntTet1 | isIntTet2;
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
        (sum(abs(vSign(face)), 2) >= 10^(-12));% & (sum(abs(vSign(face)), 2) > 0);
    c2 = ct + sum(isPolyTriFace);
    tface((ct+1):c2,:) = face(isPolyTriFace,:);
    tface2elem((ct+1):c2,:) = intIdx2cube(isPolyTriFace); % the indices of the polyhedron
    ct = c2;
    % note that only interior elements are treated
    
    % find the triangle faces on interface with normal points to exterior
    % 1. face is on the boundary and
    % 2. all vertices are on the interface
    isInterface = (sum(abs(vSign(face)), 2) <= 10^(-12));
    
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

% plot to check
%trisurf(tface(c2old+1:c2,:),node(:,1),node(:,2),node(:,3))

tface((ct+1):end,:) = [];    % remove empty meomory
tface2elem((ct+1):end) = [];
face2elemLoc = [ones(c2old,1);2*ones(c2-c2old,1)]; % face2elemLoc is the same as face location


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
meshI.vSign = vSign;

end