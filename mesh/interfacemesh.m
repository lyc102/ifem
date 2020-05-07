function [node,elem,interfaceData] = interfacemesh(box,phi,h)
%% INTERFACEMESH generated interface fitted mesh
%
%
%
%
%  Example
%     box = [ -1, 1, -1, 1];
%     h = 0.1;
%     phi = @(p) sum(p.^2, 2) - 0.5.^2;
%     [node,elem,interfaceData] = interfacemesh(box,phi,h);
%
% To Do: add option to showmesh. M-lint
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
% compute the phi value at each vertex
phiValue = phi(node);
phiValue(abs(phiValue) < 0.1*h^2) = 0;  % treat nearby nodes as on the interface
vSign = msign(phiValue);

%% Step 1: Find points near interface 
% Find the intersection points between edges and the interface
isCutEdge = (vSign(edge(:,1)).* vSign(edge(:,2))<0);
A = node(edge(isCutEdge,1),:);
B = node(edge(isCutEdge,2),:);
cutNode = findintersectbisect(phi,A,B);
Ncut = size(cutNode, 1);
vSign(N+1:N+Ncut) = 0;

% find interface elem and nodes
isInterfaceElem = false(NT,1);  
isInterfaceElem(edge2elem(isCutEdge,[1,2])) = true;
isInterfaceElem(sum(abs(vSign(elem)), 2) < 3) = true; % 1 or 2 vertices on interface
isInterfaceNode = false(N,1);
isInterfaceNode(elem(isInterfaceElem,:)) = true;

% add centers of special elements
%  0 - -1    -1 - 0     0 -  0
%  1 -  0     0 - 1     1 - -1
isSpecialElem = (sum(vSign(elem),2) == 0) & (sum(abs(vSign(elem)),2) == 2);
selem = elem(isSpecialElem,:);
% eliminate the third case
isSpecialElem = (vSign(selem(:,1)).*vSign(selem(:,3)) == -1) | ...
                (vSign(selem(:,2)).*vSign(selem(:,4)) == -1);
selem = selem(isSpecialElem,:);
auxPoint = (node(selem(:,1),:) + node(selem(:, 3),:))/2.0;
Naux = size(auxPoint,1);

% add intersected points and aux points
node = [node; cutNode; auxPoint];
interfaceNode = [node(isInterfaceNode,:); cutNode; auxPoint];
interfaceNodeIdx = [find(isInterfaceNode); N+(1:Ncut)'; N+Ncut+(1:Naux)'];
vSign(N+Ncut+1:N+Ncut+Naux) = 0;

%% Step 2: generate a Delaunay triangulation of interface points
% construct the Delaunay triangulation of interfaceNode
% different versions of matlab using different delaunay triangulation
matlabversion = version();
if str2double(matlabversion(end-5:end-2)) <= 2013
    DT = DelaunayTri(interfaceNode); %#ok<*DDELTRI>
    tElem = DT.Triangulation;
else
    DT = delaunayTriangulation(interfaceNode);
    tElem = DT.ConnectivityList;
end
tElem = fixorder(interfaceNode,tElem); % correct the orientation

%% Step 3: Post-processing
% get rid of the unnecessary triangles
NI = sum(isInterfaceNode);  % number of vertices of interface elements
isUnnecessaryElem = (sum(tElem<=NI,2) == 3); % if three vertices are vertices
tElem = tElem(~isUnnecessaryElem,:); % of interface elements, it is not needed
tElem = interfaceNodeIdx(tElem);  % map interfaceNode index to node index
% get the remainding quad elems
sElem = elem(~isInterfaceElem,:);
% merge into one triangulation
elem = [tElem; sElem(:,[2 3 1]); sElem(:,[4 1 3])];
% find interface edges and nodes on the interface
isExteriorElem = false(size(tElem,1),1);
center = (node(tElem(:,1),:) + node(tElem(:,2),:) + node(tElem(:,3),:))/3;
isExteriorElem(phi(center)>0) = true;
[~,bdEdge] = findboundary(tElem(isExteriorElem,:));
isInterfaceEdge = ((vSign(bdEdge(:,1)) == 0) & vSign(bdEdge(:,2)) == 0);
interface = bdEdge(isInterfaceEdge,:);
interfaceNode = find(vSign == 0);

%% Generate interfaceData
interfaceData.vSign = vSign;
interfaceData.tElem = tElem;
interfaceData.sElem = sElem;
interfaceData.interface = interface;
interfaceData.interfaceNode = interfaceNode;