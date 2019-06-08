%% RATE OF CONVERGENCE OF ADAPTIVE FINITE ELEMENT METHOD USING WG ELEMENT
%
% This example is to show the rate of convergence of the lowest order
% finite element approximation of the second order elliptic equation.
%
% # Lshape problem.
% # Kellogg problem.

clear variables;
close all;

%% Lshape problem
[node,elem] = squaremesh([-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
pde = Lshapedata;
format shorte
option.L0 = 1;
option.maxIt = 50;
option.maxN = 1e4;
option.printlevel = 1;
option.elemType = 'WG';
option.plotflag = 1;
err = afemPoisson(mesh,pde,option);
figure;
showrate2(err.N,err.H1,10,'k-*','||Du-Du_h||',err.N,err.eta,10,'-+','eta');

%% Kellogg problem
[node,elem] = squaremesh([-1 1 -1 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
pde = Kelloggdata;
% singular node
distance = (node(:,1)-0).^2 + (node(:,2)-0).^2;
[mindist,singularnode] = min(distance); %#ok<ASGLU>
pde.singularnode = singularnode;
option.L0 = 0;
option.maxIt = 100;
option.maxN = 1e4;
option.theta = 0.1;
option.elemType = 'WG';
option.plotflag = 1;
% option.printlevel = 1;
err = afemPoisson(mesh,pde,option);
figure;
showrate2(err.N,err.H1,10,'k-*','||Du-Du_h||',err.N,err.eta,10,'-+','eta');