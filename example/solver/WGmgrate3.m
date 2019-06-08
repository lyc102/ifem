%% WGMGRATE3 Test convergence of multigrid methods for weak Galerkin methods in 3D
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
option.printlevel = 1;
option.plotflag = 1;
option.rateflag = 0;
option.dispflag = 0;
colname = {'#Dof','Steps','Time'};

%% Example: cube mesh and Poisson equation
pde = sincosdata3;
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);
option.maxN = 3e6;
showmesh3(node,elem,[],'Facecolor','w','FaceAlpha',0.4);
bdFlag = setboundary3(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

option.reducesystem = 1; % Solve the reduced system
[err,time,solver] = femPoisson3(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');

option.reducesystem = 0; % Solve the original system
[err,time,solver] = femPoisson3(mesh,pde,option);
disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');

%% Example: cube mesh and jump coefficients
pde = jumpmgdata2;
global epsilon
option.maxN = 3e6;
for k = -4:2:4
    epsilon = 10^k;
    fprintf('epsilon = %8.2e \n',epsilon);
    [node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);
    bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | (x==-1)');
    mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
    option.reducesystem = 1; % Solve the reduced system
    [err,time,solver] = femPoisson3(mesh,pde,option);
    disptable(colname,solver.N,[],solver.itStep,[],solver.time,'%4.2g');
end