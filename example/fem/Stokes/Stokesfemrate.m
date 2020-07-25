%% FINITE ELEMENTS FOR STOKES EQUATIONS IN 2D
%
% This example is to show the convergence of various finite element
% approximation of the Stokes equation on the unit square:
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%                              u = g_D  on \Gamma.
%
% with the pure Dirichlet boundary condition. The solver is based on a DGS
% type smoother. 
%
% Reference
%
% L. Chen. Finite element methods for Stokes equations. Course notes.

clear variables; 
close all;

%% Setting
% mesh
[node,elem] = squaremesh([0,1,0,1],0.25);
showmesh(node,elem);
[node,elem] = uniformrefine(node,elem);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
figure; showmesh(node,elem); pause(0.5);
% pde
pde = Stokesdata1; 
% options
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.solver = 'mg';

%% CR-P0 element
disp('CR-P0')
option.elemType = 'CRP0';
femStokes(mesh,pde,option);

%% P2-P0 element
disp('P2-P0')
option.elemType = 'P2P0';
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
option.solver = 'asmg';
femStokes(mesh,pde,option);

%% P2-P1 element
disp('P2-P1')
option.elemType = 'P2P1';
femStokes(mesh,pde,option);
