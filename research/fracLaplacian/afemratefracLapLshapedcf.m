close all; clear all;

%% Setting of the problem
global s
pde = fracLapdata8; % f = 1 for |x - (-0.5,0.5)| < 0 .25
pde.L = 1;
option.theta = 0.3;
option.estType = 'star';
option.maxN = 2e5;
option.solver = 'mg';
option.gNquadorder = 5;
option.tol = 1e-6;
[node,elem] = squaremesh([-1,1,-1,1],0.5);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');

%% s = 0.2
s = 0.2;
option.maxIt = 17;
afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.4
s = 0.4;
option.maxIt = 19;
afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.6
s = 0.6;
option.maxIt = 19;
afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.8
s = 0.8;
option.maxIt = 19;
afemfracLap(node,elem,pde,bdFlag,option);