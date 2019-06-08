close all; clear all;

%% Setting of the problem
global s
pde = fracLapdata7; % f = 1 for |x|<1/2
pde.L = 1;
option.theta = 0.3;
option.estType = 'star';
option.maxIt = 16;
option.maxN = 3e4;
option.solver = 'mg';
option.gNquadorder = 5;
option.tol = 1e-6;
[node,elem] = squaremesh([-1 1 -1 1],0.5);
bdFlag = setboundary(node,elem,'Dirichlet');

%% s = 0.2
s = 0.2;
afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.4
s = 0.4;
afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.6
s = 0.6;
afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.8
s = 0.8;
afemfracLap(node,elem,pde,bdFlag,option);