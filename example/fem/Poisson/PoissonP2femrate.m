%% RATE OF CONVERGENCE OF QUADRATIC ELEMENT FOR POISSON EQUATION
%
% This example is to show the rate of convergence of quadratic finite
% element approximation of the Poisson equation on the unit square with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
% - Robin boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem PoissonP2femrate">PoissonP2femrate</a>.
%
% See also Poissonfemrate, Poissonafemrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
mesh = struct('node',node,'elem',elem);
option.L0 = 2;
option.maxIt = 4;
option.printlevel = 1;
option.plotflag = 1;
option.elemType = 'P2';

%% Non-empty Dirichlet boundary condition.
pde = sincosdata;
mesh.bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
femPoisson(mesh,pde,option);

%% Pure Neumann boundary condition.
option.plotflag = 0;
pde = sincosNeumanndata;
mesh.bdFlag = setboundary(node,elem,'Neumann');
femPoisson(mesh,pde,option);

%% Pure Robin boundary condition.
option.plotflag = 0;
pde = sincosRobindata;
mesh.bdFlag = setboundary(node,elem,'Robin');
femPoisson(mesh,pde,option);

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (2nd order) and L2-norm
% (3rd order) is observed. The order of ||DuI-Duh|| is 2.5 order and
% thus superconvergence exists.
%
% MGCG converges uniformly in all cases.