% To Do: add more examples

clear all;
close all;

box = [ -1, 1, -1, 1];
hmax = 0.1;
phi = @(p) sum(p.^2, 2) - 0.5.^2;
[node,tElem,sElem,bdEdge, interfaceNodeIdx] = interfacemesh(phi,box,hmax);