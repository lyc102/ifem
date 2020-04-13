function  [node,face,face2elem,interfaceData] = interfacemesh3(cube,surface,h)
%% INTERFACEMESH3 generates the mesh structure for elliptic interface problems
%
%   [node,face,face2elem,interfaceData] = INTERFACEADAPTIVEMESH(cube, phi, h), where 
%    - cube = [xmin, xmax, ymin, ymax, zmin, zmax],  represent an cube domain 
%       enclosing the interface; 
%    - surface class 
%      - surface.phi is a level set function defining the interface as follows
%       {p | phi(p)=0} represents the interface
%       {p | phi(p) < 0} represents the interior domain enclosed by the interface
%       {p | phi(p) > 0} rpepresents the exterior domain of the interface. 
%      - surface.number number of different level set functions 
%      - surface.phi_i  i-th level set function 
%    - h is the mesh size of the Cartesian grid on cube.
%
%   Example: One sphere
%       cube = [-1 1 -1 1 -1 1];
%       surface.phi = @(p) sum(p.^2, 2) - 0.5.^2;
%       h = 0.125;
%       [node,face,face2elem,interfaceData] = interfacemesh3(cube,surface,h);
%       showmesh(node,interfaceData.interface);
%
%   Example: Two spheres
%       cube = [-1 1 -1 1 -1 1]*2;
%       surface = twospherephi;
%       h = 0.125;
%       [node,face,face2elem,interfaceData] = interfacemesh3(cube,surface,h);
%       showmesh(node, interfaceData.interface);
%
%   Author: Huayi Wei <weihuayi@xtu.edu.cn>, based on discussion with Long Chen.
%
%   See also: interfacemesh, interfaceadaptivemesh
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Generate the Cartesian mesh on a cube domain with size h
[node,elem] = cubehexmesh(cube, h);
N = size(node, 1); % number of nodes
NC = size(elem,1); % number of small cubes

%% Evaluate sign of vertices and find interface elements
% find the sign function of all vertices
phiValue = surface.phi(node);
vSign = msign(phiValue);
e2vSign = vSign(elem);     % element to vertices sign

% find cut element
eta = max(e2vSign,[],2).*min(e2vSign,[],2);
isCutElem  = (eta < 0 | sum(abs(e2vSign),2) <= 5); 
% eta < 0: elements cross interface
% sum(abs(e2vSign),2) <= 5: at least 3 vertices on the interface
cutElem = elem(isCutElem,:);
% find more cut element
c = (node(elem(:,1),:)+node(elem(:,8),:))/2;
dist = surface.phi(c);
elemSign = sum(e2vSign,2);
isSameElem = (abs(elemSign) == 8);
isMoreCutElem = (abs(dist)< h/2 - eps) & isSameElem;
% distance of center to surface < h/2 but all vertices have the same sign
moreCutElem = elem(isMoreCutElem,:);
isallCutElem = isCutElem | isMoreCutElem;
allCutElem = elem(isallCutElem,:);
 
% index of element
elemIdx = zeros(NC,1,'int8');
elemIdx((eta >= 0) & (max(e2vSign,[],2) == 1)) = 1;   % 1: exterior element
elemIdx((eta >= 0) & (min(e2vSign,[],2) == -1)) = -1; % -1: interior element
elemIdx(isallCutElem) = 0; % 0: interface element

%% Find intersection points
% find the cut edge
localEdge = [1 5; 2 6; 4 8; 3 7; 1 2; 3 4; 5 6; 7 8; 1 4; 2 3; 5 8; 6 7];
edgeMatrix = sparse(cutElem(:,localEdge), cutElem(:,localEdge(:,[2,1])), 1, N, N);
[I,J,~] = find(triu(edgeMatrix));
edge = [I, J];
isCutEdge = vSign(edge(:,1)).*vSign(edge(:, 2)) <  0; 
cutEdge = edge(isCutEdge, :);

% find the cut points
A = node(cutEdge(:,1),:);
B = node(cutEdge(:,2),:);
cutPoint = findintersectbisect(@(p)surface.phi(p),A,B);
nCutPt = size(cutPoint,1); 
N1 = N + nCutPt;
node(N+1:N1,:) = cutPoint; % add the coordinates of the intersection points
vSign(N+1:N1) = 0;

% add more intersection points in elements with the same sign
localEdge = [1 5; 2 6; 4 8; 3 7; 1 2; 3 4; 5 6; 7 8; 1 4; 2 3; 5 8; 6 7];
edgeMatrix = sparse(moreCutElem(:,localEdge), moreCutElem(:,localEdge(:, [2,1])), 1, N, N);
[I,J,~] = find(triu(edgeMatrix));
moreEdge = [I, J];
moreCutNode = zeros(size(moreEdge,1),3); % preallocate
if ~isfield(surface,'number')
    surface.number = 0;
end
nMore = 1;
for i = 1:surface.number % compute intersection for each level set
    isSpecialCutEdge = surface.phi_i(node(moreEdge(:, 1), :),i)...
                     .*surface.phi_i(node(moreEdge(:, 2), :),i) < 0;
    cutEdge = moreEdge(isSpecialCutEdge, :);
    A = node(cutEdge(:,1),:);
    B = node(cutEdge(:,2),:);
    cutNode = findintersectbisect(@(p)surface.phi_i(p, i), A, B);
    % add only nodes on the interface
    isSurfaceNode = abs(surface.phi(cutNode)) < 1e-8;
    cutNode = cutNode(isSurfaceNode, :);
    if ~isempty(cutNode)
        nCutNode = size(cutNode,1); 
        moreCutNode(nMore:nMore + nCutNode-1,:) = cutNode;
        nMore = nMore + nCutNode-1;
    end
end
if nMore > 1  % more cute points are found
    node(N1+1:N1+nMore, :) = moreCutNode(1:nMore,:);
    vSign(N1+1:N1+nMore) = 0;
end

%% Find the special square face with diagonal vertex on the interface
T = auxstructurehex(allCutElem);
face = T.face;
face2elem = double(T.face2elem);
isSpecialFace = (face2elem(:, 1) ~= face2elem(:, 2)) & ...
                ((sum(vSign(face),2) == 0 & sum(abs(vSign(face)),2) == 2));
% add the center into the interface points
auxPoint = (node(face(isSpecialFace,1),:)+node(face(isSpecialFace,2),:)...
          + node(face(isSpecialFace,3),:)+node(face(isSpecialFace,4),:))/4.0;
N2 = N1 + nMore+ size(auxPoint,1);
node(N1+nMore+1:N2, :) = auxPoint; % add the coordinates of the aux points
vSign(N1+nMore+1:N2, :) = 0;

%% Find the boundary faces of all cut elements which are squares.
sface = double(T.bdFace);           % node idx is global
sface2elem = double(T.bdFace2elem); % elem idx is to allCutElem
localidx2globalidx = find(isallCutElem);
sface2elem = localidx2globalidx(sface2elem); % map to idx of elem
eta = max(vSign(sface), [] ,2) == 1;  % outerior faces
sface2elem(eta) = sface2elem(eta) + NC; % outerior polyhedron idx
clear T;

%% Construct Delaunay triangulation
isInterfaceNode = false(size(node,1), 1);
isInterfaceNode(allCutElem(:)) = true; % include the vertices of allCutElem 
isInterfaceNode(N+1:end) = true;    % and the cut points and the aux points
interfaceNode = node(isInterfaceNode,:);
% different versions of matlab using different delaunay triangulation
matlabversion = version();
if str2double(matlabversion(end-5:end-2)) <= 2013
    DT = DelaunayTri(interfaceNode); %#ok<*DDELTRI>
    tetElem = DT.Triangulation;
else
    DT = delaunayTriangulation(interfaceNode);
    tetElem = DT.ConnectivityList;
end
tetElem = fixorder3(interfaceNode, tetElem);
localidx2globalidx = find(isInterfaceNode); % tetElem points to interfaceNode
tetElem = localidx2globalidx(tetElem);  % map to the global index of node
bc = (node(tetElem(:, 1), :) + node(tetElem(:, 2), :) ...
    + node(tetElem(:, 3), :) + node(tetElem(:, 4), :))/4.0;
idx2cube = getCubeIdx(bc, cube, h); % Get the cube index of each tetElem
tetElem = tetElem(isallCutElem(idx2cube), :); % Just keep the tets in cut elem
idx2cube = idx2cube(isallCutElem(idx2cube));  % keep cube idx in cut elements

% There might exist some bad tet elements whose vertices are in
% different cubes. We just need to keep the tet in every cut elem. So
% here we get rid of them. 
NT = size(tetElem,1);
X = reshape(node(tetElem,1), NT, 4);
Y = reshape(node(tetElem,2), NT, 4);
Z = reshape(node(tetElem,3), NT, 4);
isBadTet =  (max(X,[],2) - min(X,[],2)) > 2*h - eps ... % d > 2h means
          | (max(Y,[],2) - min(Y,[],2)) > 2*h - eps ... % in differnt cube
          | (max(Z,[],2) - min(Z,[],2)) > 2*h - eps;
tetElem = tetElem(~isBadTet,:);
idx2cube = idx2cube(~isBadTet);

%% Get triangular faces for interrior elements and interface
localFace = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
NT = size(tetElem,1);
tface = zeros(4*NT, 3);
tface2elem = zeros(4*NT, 1);
interface = zeros(4*NT, 3);

% find the interior tet elements
isIntTet = min(vSign(tetElem),[], 2) == -1;
intTet = tetElem(isIntTet,:);
% find the corresponding cube indices
intIdx2cube = idx2cube(isIntTet); 
% find triangular interface
T = auxstructure3(intTet);
neighbor = T.neighbor;
clear T;
tmp = (1:size(intTet, 1))';
ct = 0;
ci = 0;
for i = 1:4
    face = intTet(:, localFace(i,:));
    % find the triangle faces of polyhedron
    % 1. face and its neighbor associated to different cubes, and
    % 2. face is not on the boundary of all cut elements (which are squares)
    isPolyTriFace = (neighbor(:, i) ~= tmp) & (intIdx2cube ~= intIdx2cube(neighbor(:,i)));
    c2 = ct + sum(isPolyTriFace);
    tface((ct+1):c2,:) = face(isPolyTriFace,:);
    tface2elem((ct+1):c2,:) = intIdx2cube(isPolyTriFace); % the indices of the polyhedron
    % note that only interior elements are treated
    
    % find the triangle faces on interface with normal points to exterior
    % 1. face is on the boundary and
    % 2. all vertices are on the interface
    isInterface = (neighbor(:, i) == tmp) & (sum(abs(vSign(face)), 2) == 0);
    ct = c2;
    c2 = ct + sum(isInterface);
    % add to polyhedron
    tface((ct+1):c2,:) = face(isInterface,:); % the interface tri faces.
    tface2elem((ct+1):c2) = intIdx2cube(isInterface); % the indices of the polyhedron
    ct = c2;
    % add to interface    
    c4 = ci + sum(isInterface);
    interface((ci+1):c4,:) = face(isInterface,:);
    ci = c4;
end
interface((ci+1):end,:) = [];

% find the interface 
center = (node(intTet(:,1),:)+node(intTet(:,2),:)+node(intTet(:,3),:)+node(intTet(:,4),:))/4;
isIntElem = surface.phi(center)<0;
[~,bdface] = findboundary3(intTet(isIntElem,:));

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
    isPolyTriFace = (neighbor(:, i) ~= tmp) & (extIdx2cube ~= extIdx2cube(neighbor(:,i)));
    c2 = ct + sum(isPolyTriFace);
    tface((ct+1):c2,:) = face(isPolyTriFace,:);
    tface2elem((ct+1):c2,:) = extIdx2cube(isPolyTriFace) + NC; 
    % index of exterior polyhedron is append to the end of elem
    ct = c2;
    % find the triangle faces on interface
    % 1. face is on the boundary and
    % 2. all vertices are on the interface
    isInterface = (neighbor(:, i) == tmp) & (sum(abs(vSign(face)), 2) == 0);    
    c2 = ct + sum(isInterface);
    tface((ct+1):c2,:) = face(isInterface,:);
    tface2elem((ct+1):c2) = extIdx2cube(isInterface) + NC;
    ct = c2;
end
% check connectedness of polyhedron

% generate face and face2elem
tface((ct+1):end,:) = [];    % remove empty meomory
tface2elem((ct+1):end) = [];
face2elem = [tface2elem; sface2elem];
face = [tface,zeros(size(tface,1),1); sface];

%% Generate interfaceData
interfaceData.elemIdx = elemIdx;
interfaceData.vSign = vSign;
interfaceData.interface = interface;
interfaceData.cutElem = allCutElem;