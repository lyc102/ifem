%% CONVERGENCE OF FAST SOLVERS FOR MIXED FINITE ELEMENT METHOD (RT0-P0) FOR DARCY'S EQUATIONS
%
% This example is to show the convergence of fast solvers of mixed finite
% element (RT0-P0) approximation of the Darcy's equations.
%
% Reference 
%
% L. Chen. Multigrid Methods for Constrained Minimization Problems and
% Application to Saddle Point Problems. Submitted,2014.

close all
clear variables
%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
mesh = struct('node',node,'elem',elem);
option.L0 = 1;
option.maxIt = 4;
option.printlevel = 1;
option.elemType = 'RT0';
option.refType = 'red';
% option.solver = 'uzawapcg';
% option.solver = 'tri';
option.solver = 'mg';

%% Poisson
pde = Darcydata0;
disp('Poisson: uniform grid')
mesh.bdFlag = setboundary(node,elem,'Neumann');
mfemDarcy(mesh,pde,option);

%% Anisotropic tensor
disp('Anisotropic tensor: uniform grid')
[node,elem] = squaremesh([0,1,0,1],0.25); 
pde = Darcydata4;
option.dquadorder = 3;
mesh.bdFlag = setboundary(node,elem,'Neumann');
mfemDarcy(mesh,pde,option);

%% Jump tensor
disp('Jump tensor: uniform grid')
pde = Darcydata3;
mesh.bdFlag = setboundary(node,elem,'Neumann');
mfemDarcy(mesh,pde,option);

%% Jump tensor with distortion of grids
disp('Jump tensor: distorted grid')
% perturbe the grid
[bdNode,bdEdge,isBdNode] = findboundary(elem);
isIntNode = ~isBdNode;
mesh.node(isIntNode,:) = mesh.node(isIntNode,:) + 0.25*0.4*rand(sum(isIntNode),2);
pde = Darcydata3;
mesh.bdFlag = setboundary(node,elem,'Neumann');
mfemDarcy(mesh,pde,option);
