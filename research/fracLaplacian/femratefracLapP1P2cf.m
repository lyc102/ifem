%% Setting of the problem
global s
pde = fracLapdata9;
% pde = checkboarddata;
pde.L = 1;
option.maxIt = 4;
option.maxN = 1e6;
option.elemType = 'P1P2';
option.solver = 'mg';
option.gNquadorder = 4;
option.tol = 1e-6;
[node,elem] = squaremesh([-1 1 -1 1],0.5);
bdFlag = setboundary(node,elem,'Dirichlet');

%% s = 0.2
s = 0.2;
femfracLap(node,elem,pde,bdFlag,option);

%% s = 0.4
s = 0.4;
femfracLap(node,elem,pde,bdFlag,option);

%% s = 0.6
s = 0.6;
femfracLap(node,elem,pde,bdFlag,option);

%% s = 0.8
s = 0.8;
femfracLap(node,elem,pde,bdFlag,option);