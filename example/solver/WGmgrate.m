%% WGMGRATE Test convergence of multigrid methods for weak Galerkin method
%
% Reference
%
% An auxiliary space multigrid preconditioner for the weak Galerkin method.
% By Long Chen, Junping Wang, Yanqiu Wang, and Xiu Ye. Computers and
% Mathematics with Applications, 70(4):330?344 2015.
% doi:10.1016/j.camwa.2015.04.016


%% Options
clear variables; 
close all;
option.maxIt = 4;
option.elemType = 'WG';
option.solver = 'mg';
option.smoothingstep = 2;
option.plotflag = 1;
option.rateflag = 0;
option.dispflag = 1;
colname = {'#Dof','Steps','Time'};

%% Example: circle mesh and Poisson equation
pde = sincosdata;
[node,elem] = circlemesh(0,0,1,0.25);
% showmesh(node,elem,'Facecolor','w');
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

option.L0 = 1;
option.maxN = 3e5;
option.reducesystem = 0; % Solve the original system
[err,time,solver] = femPoisson(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');

option.reducesystem = 1; % Solve the reduced system
[err,time,solver] = femPoisson(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');

%% Example: square mesh and variable Poisson equation
[node,elem] = squaremesh([0,1,0,1],0.25); 
% showmesh(node,elem,'Facecolor','w');
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

option.L0 = 2;
option.dquadorder = 4;
pde = oscdiffdata;

option.reducesystem = 0; % Solve the original system
[err,time,solver] = femPoisson(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');

option.reducesystem = 1; % Solve the reduced system
[err,time,solver] = femPoisson(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');

%% Example: adaptive mesh and Poisson equation with less regularity
[node,elem] = squaremesh([-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

pde = Lshapedata;
option.maxIt = 30;
option.maxN = 1e4;
option.printlevel = 0;
option.rateflag = 1;

option.reducesystem = 0; % Solve the original system
[err,time,solver] = afemPoisson(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');

option.reducesystem = 1; % Solve the reduced system
[err,time,solver,eqn,node,elem] = afemPoisson(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');
% figure; showmesh(node,elem,'Facecolor','w');