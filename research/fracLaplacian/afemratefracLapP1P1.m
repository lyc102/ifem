close all; clear all;

%% Setting of the problem
global s
pde = fracLapdata2; % smooth solution u =  sin(2*pi*x1)*sin(2*pi*x2);
pde.L = 1;
option.theta = 0.3;
% option.estType = 'star';
option.estType = 'recovery';
option.maxIt = 25;
option.maxN = 2e5;
option.solver = 'mg';
option.gNquadorder = 4;
option.tol = 1e-6;
[node,elem] = squaremesh([0 1 0 1],0.25);
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