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
% Reference
%
% M. Wang and L. Chen. Multigrid Methods for the Stokes equations
% using Distributive Gauss-Seidel Relaxations based on the Least Squares
% Commutator. Journal of Scientific Computing. 56(2): 409-431, 2013.

clear variables; 
close all;

%% Setting
[node,elem] = squaremesh([0,1,0,1],0.125);
% [node,elem] = circlemesh(0,0,1,0.25);
% [node,elem] = uniformrefine(node,elem);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
figure;
showmesh(node,elem); pause(0.5);
% pde
pde = Stokesdata1; 
%  pde = StokesZulehnerdata;
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.plotflag = 0;
option.rateflag = 0;

%% MG options
option.solver = 'asmg';
option.smoothingstep = 2;
option.smootherbarSp = 'SGS';

%% CR-P0 element
disp('CR-P0')
option.elemType = 'CRP0';
femStokes(mesh,pde,option);

%% P2-P0 element
disp('P2-P0')
option.elemType = 'P2P0';
femStokes(mesh,pde,option);

%% P2-P1 element
disp('P2-P1')
option.elemType = 'P2P1';
option.solver = 'asmg';
femStokes(mesh,pde,option);

%% isoP2-P0 element
disp('isoP2-P0')
option.elemType = 'isoP2P0';
femStokes(mesh,pde,option);

%% isoP2-P1 element
disp('isoP2-P1')
option.elemType = 'isoP2P1';
femStokes(mesh,pde,option);

%% P1b-P1 element
disp('P1b-P1')
option.elemType = 'P1bP1';
femStokes(mesh,pde,option);