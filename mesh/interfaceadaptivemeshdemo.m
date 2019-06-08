clear all; close all;
% initial mesh
node = [-1 -1; 1 -1; 1 1; -1 1];
elem = [2 3 1; 4 1 3];
% interface mesh
% phi is the level set function
[node,elem,phiValue] = interfaceadaptivemesh(node,elem,@(p)phi1(p,12));
showmesh(node,elem);
idx = (abs(phiValue) < eps);
findnode(node,idx,'noindex');
