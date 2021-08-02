function [node,elem,interface] = interfacemeshdoc(box,phi,h)
%% INTERFACEMESHDOC doc for interfacemesh
%
% [node,elem,interfaceData] = INTERFACEMESH(box,phi,h) generates an
% interface-fitted mesh [node,elem]. The box is a rectangle containing the
% interface and the interface is defined as the zero level set of function
% phi. The parameter h is the size of the backgroud uniform mesh.
%
% Information on the discrete interface is stored in the output
% interface which includes:
%   - interface.vSign: sign of phi at vertices;
%        1: outside;  -1: inside;  0: on the interface
%   - interface.tElem: traingles near the interface;
%   - interface.sElem: traingles away from the interface;
%   - interface.edge: edges approximating the interface;
%   - interface.node: vertices on the interface;
%
%  Examples
%   % circle
%     box = [ -1, 1, -1, 1];
%     h = 0.1;
%     phi = @(p) sum(p.^2, 2) - 0.5.^2;
%     [node,elem,interfaceData] = interfacemeshdoc(box,phi,h);
%
%    % heart
%     box = [ -1, 1, -1, 1];
%    [node,elem,interface] = interfacemeshdoc(box,@phiheart,0.05);
%    showmesh(node,elem);
%    findedge(node,interface.edge,'all','noindex','draw');
%
%
% Author: Huayi Wei <weihuayi@xtu.edu.cn>, and Long Chen.
%
%   Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Construct the initial structure mesh
[node,elem,T] = squarequadmesh(box,h);
edge = T.edge;
edge2elem = T.edge2elem;
clear T;
N = size(node,1);
NT = size(elem,1);
showmesh(node,elem);
pause;

%% Step 1: Find points near interface 
% compute the phi value at each vertex
phiValue = phi(node);
phiValue(abs(phiValue) < eps*h) = 0;  % treat nearby nodes as on the interface
vSign = sign(phiValue);
findnode(node,vSign==0);

% Find the intersection points between edges and the interface
isCutEdge = (vSign(edge(:,1)).* vSign(edge(:,2))<0);
findedge(node,edge,isCutEdge,'noindex','draw');
A = node(edge(isCutEdge,1),:);
B = node(edge(isCutEdge,2),:);
cutNode = findintersectbisect(phi,A,B);
Ncut = size(cutNode, 1);
vSign(N+1:N+Ncut) = 0;
findnode(cutNode,'all','noindex','color','k','MarkerSize',12);
pause;

% find interface elem and nodes
isInterfaceElem = false(NT,1);  
isInterfaceElem(edge2elem(isCutEdge,[1,2])) = true;
isInterfaceElem(sum(abs(vSign(elem)), 2) < 3) = true; % 2 vertices on interface
findelem(node,elem,isInterfaceElem,'noindex');        

isInterfaceNode = false(N,1);
isInterfaceNode(elem(isInterfaceElem,:)) = true;
findnode(node,isInterfaceNode,'noindex');
pause;

% add centers of special elements (two vertices on the interface)
%  0 - -1    -1 - 0     0 -  0      1 - 0
%  1 -  0     0 - 1     1 - -1      1 - 0
isSpecialElem = (sum(abs(vSign(elem)),2) <= 2); % two or more vertices on interfaces
% isSpecialElem = (sum(vSign(elem),2) == 0) & (sum(abs(vSign(elem)),2) == 2);
specialElem = elem(isSpecialElem,:);
auxPoint = (node(specialElem(:,1),:) + node(specialElem(:, 3),:))/2.0;
Naux = size(auxPoint,1);
vSign(N+Ncut+1:N+Ncut+Naux) = sign(phi(auxPoint));
% find the first two cases
isDiagInterface = (vSign(specialElem(:,1)).*vSign(specialElem(:,3)) == -1) | ...
                  (vSign(specialElem(:,2)).*vSign(specialElem(:,4)) == -1);
vSign(N+Ncut+find(isDiagInterface)) = 0;            
findnode(auxPoint,'all','noindex','color','m');
pause;

% add cut points and aux points
node = [node; cutNode; auxPoint];
nearInterfaceNode = [node(isInterfaceNode,:); cutNode; auxPoint];
interfaceNodeIdx = [find(isInterfaceNode); N+(1:Ncut)'; N+Ncut+(1:Naux)'];
% figure(2);
% showmesh(node,elem);
findnode(nearInterfaceNode,'all','noindex');
pause;

%% Step 2: generate a Delaunay triangulation of interface points
% construct the Delaunay triangulation of interfaceNode
% different versions of matlab using different delaunay triangulation
matlabversion = version();
if str2double(matlabversion(end-5:end-2)) <= 2013
    DT = DelaunayTri(nearInterfaceNode); %#ok<*DDELTRI>
    tElem = DT.Triangulation;
else
    DT = delaunayTriangulation(nearInterfaceNode);
    tElem = DT.ConnectivityList;
end
tElem = fixorder(nearInterfaceNode,tElem); % correct the orientation
% showmesh(node,tElem);

%% Step 3: Post-processing
% get rid of the unnecessary triangles
NI = sum(isInterfaceNode);  % number of pts of the background mesh near interface
haveNewPts = (sum(tElem > NI,2) > 0); % triangles containing new added vertices
% vSignNearInterface = zeros(size(nearInterfaceNode,1),1);
% vSignNearInterface(1:sum(isInterfaceNode)) = vSign(isInterfaceNode);
% haveInterfacePts = (sum(abs(vSignNearInterface(tElem)),2) < 3);
tElem = tElem(haveNewPts,:); 
% tElem = tElem(haveInterfacePts,:); 
tElem = interfaceNodeIdx(tElem);  % map interfaceNode index to node index
showmesh(node,tElem,'Facecolor','y');
findnode(node,vSign==0,'noindex');
pause;

% get the remainding quad elems
sElem = elem(~isInterfaceElem,:);
% merge into one triangulation
elem = [tElem; sElem(:,[2 3 1]); sElem(:,[4 1 3])];
T = auxstructure(tElem);
showmesh(node,elem); hold on;
showmesh(node,tElem,'Facecolor','y');
isInterfaceEdge = ((vSign(T.edge(:,1)) == 0) & vSign(T.edge(:,2)) == 0);
interfaceEdge = T.edge(isInterfaceEdge,:);
nearInterfaceNode = find(vSign == 0);
findedge(node,interfaceEdge,'all','noindex','draw');
pause;

%% Generate interfaceData
interface.vSign = vSign;
interface.tElem = tElem;
interface.sElem = sElem;
interface.edge = interfaceEdge;
interface.node = nearInterfaceNode;