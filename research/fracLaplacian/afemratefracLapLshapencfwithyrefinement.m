close all; clear all;

%% Setting of the problem
global s
pde = fonedata; % f = 1;
pde.L = 1;
option.theta = 0.3;
option.estType = 'star';
option.maxIt = 20;
option.maxN = 2e5;
option.solver = 'mg';
option.yrefinement = 10; % ratio of h_y/h_x
option.gNquadorder = 5;
option.tol = 1e-6;
[node,elem] = squaremesh([-1,1,-1,1],0.5);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
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