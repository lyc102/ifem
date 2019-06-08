%% MULTIGRID OF FOR THE STOKES EQNS IN 2D
%
% This example is to show the convergence of multigrid methods for various
% finite element approximation of the Stokes equation on the unit square:
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%                              u = g_D  on \Gamma.
%
% with the pure Dirichlet boundary condition.
%

%% Setting
clear variables; 
close all;
% [node,elem] = squaremesh([0,1,0,1],0.25);
[node,elem] = circlemesh(0,0,1,0.25); 
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% option.mesh = 'circle';
% load Lshapemesh
% load Lshapeunstructure
% load flowpastcylindermesh
showmesh(node,elem);
% pde = Stokesdata1; 
pde = StokesZulehnerdata;
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.plotflag = 0;
option.rateflag = 0;

%% MG options
option.solver = 'mg';
% option.printlevel = 2;
option.smoothingstep = 2;
% option.smootherbarSp = 'SGS';

%% RT0-P0
disp('RT0-P0')
option.elemType = 'RT0';
% option.refType = 'bisect';
femStokesHdiv(mesh,pde,option);

%% BDM1B-P0
disp('BDM1B-P0')
option.elemType = 'BDM1';
femStokesHdiv(mesh,pde,option);