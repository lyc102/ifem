---
permalink: /mesh/interfacemeshdoc/
title: "Interface Mesh Generator"
sidebar:
    nav: mesh
---



# Interface-fitted Mesh Generator: Two Dimensions

We explain a simple mesh generator for generating an interface-fitted mesh in two dimensions. 

**Algorithm:  2D Interface-fitted Mesh Generation**

0. Generate a uniform Cartesian mesh of the box

1. Find points near the interface. 
    - cut points
    - mesh points near or on the interface
    - centers of special elements.

2. Construct a Delaunay triangulation on these points.

3. Post processing. 
    - remove triangles away from the interface
    - merge all uncut elements
    
    

**Reference**: L. Chen, H. Wei, M. Wen. [An Interface-Fitted Mesh Generator and Virtual Element Methods for Elliptic Interface Problems](https://www.math.uci.edu/~chenlong/interface_vem.html). *Journal of Computational Physics*. 334(1), 2017, 327–348.

## Generate a uniform Cartesian mesh of the box


```matlab
box = [ -1, 1, -1, 1];
h = 0.1;
[node,elem,T] = squarequadmesh(box,h);
edge = T.edge;
edge2elem = T.edge2elem;
clear T;
N = size(node,1);
NT = size(elem,1);
showmesh(node,elem);
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_2_0.png)
​    


## Find points near the interface


```matlab
phi = @phiheart;
showmesh(node,elem);
hold on;

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
findnode(cutNode,'all','noindex','color','b','MarkerSize',12);
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_4_0.png)
​    



```matlab
showmesh(node,elem);
hold on;
% find interface elem and nodes
isInterfaceElem = false(NT,1);  
isInterfaceElem(edge2elem(isCutEdge,[1,2])) = true;
isInterfaceElem(sum(abs(vSign(elem)), 2) < 3) = true; % 2 vertices on interface
findelem(node,elem,isInterfaceElem,'noindex');        

isInterfaceNode = false(N,1);
isInterfaceNode(elem(isInterfaceElem,:)) = true;
findnode(node,isInterfaceNode,'noindex','MarkerSize',12);
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_5_0.png)
​    



```matlab
showmesh(node,elem);
hold on;
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
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_6_0.png)
​    



```matlab
% add cut points and aux points
node = [node; cutNode; auxPoint];
nearInterfaceNode = [node(isInterfaceNode,:); cutNode; auxPoint];
interfaceNodeIdx = [find(isInterfaceNode); N+(1:Ncut)'; N+Ncut+(1:Naux)'];
showmesh(node,elem);
findnode(nearInterfaceNode,'all','noindex','color','r','MarkerSize',11);
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_7_0.png)
​    


## Construct a Delaunay triangulation on these points


```matlab
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
showmesh(node,tElem);
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_9_0.png)
​    


The Delaunay triangulation is a mesh for the coonvex hull of the given points. Without postprocess, it is messy. 

## Post-processing


```matlab
% get rid of the unnecessary triangles
NI = sum(isInterfaceNode);  % number of pts of the background mesh near interface
haveNewPts = (sum(tElem > NI,2) > 0); % triangles containing new added vertices
tElem = tElem(haveNewPts,:); 
tElem = interfaceNodeIdx(tElem);  % map interfaceNode index to node index
showmesh(node,tElem,'Facecolor','y');
findnode(node,vSign==0,'noindex');
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_11_0.png)
​    



```matlab
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
```


​    
![png](interfacemeshdoc_files/interfacemeshdoc_12_0.png)
​    



```matlab
%% Generate interface data
interface.vSign = vSign;
interface.tElem = tElem;
interface.sElem = sElem;
interface.edge = interfaceEdge;
interface.node = nearInterfaceNode;
```

## Further Reading

Run the following example in MATLAB and read `interfacemeshdoc` to learn more.


```matlab
box = [ -1, 1, -1, 1];
[node,elem,interface] = interfacemeshdoc(box,@phiheart,0.05);
showmesh(node,elem);
findedge(node,interface.edge,'all','noindex','draw');
```
