%% RATE OF CONVERGENCE OF CUBIC FINITE ELEMENT METHOD
%
% This example is to show the rate of convergence of cubic finite element
% approximation of the Poisson equation on the regular polygons:
%
% $$- \Delta u = f \; \hbox{in } Omega$$
%
% for the Dirichlet boundary condition.

clear variables
%% Set up problem
% PDE and Boundary condition.
pde = simpledata; % f = 1, g_D = 0
% FEM
option.elemType = 'P3';

%% Options
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.plotflag = 1;

%% Case 1: Triangle
[node,elem] = regpolygon(3,0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
showmesh(node,elem);
femPoisson(mesh,pde,option);

%% Case 2: Square
[node,elem] = squaremesh([0,1,0,1],0.25); 
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
showmesh(node,elem);
femPoisson(mesh,pde,option);

%% Case 3: Pentagon
[node,elem] = regpolygon(5,0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
showmesh(node,elem);
femPoisson(mesh,pde,option);

%% Case 4: Hexagon
[node,elem] = regpolygon(6,0.5);
showmesh(node,elem);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
showmesh(node,elem);
femPoisson(mesh,pde,option);