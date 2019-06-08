%% RATE OF CONVERGENCE OF BILINEAR FINITE ELEMENT METHOD
%
% This example is to show the rate of convergence of bilinear finite
% element approximation of the Poisson equation on the unit square with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
% - Robin boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem PoissonQ1femrate">PoissonQ1femrate</a>.
%
% See also Poissonfemrate, Poissonafemrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear variable
%% Setting
[node,elem] = squarequadmesh([0,1,0,1],1/2^3); 
mesh  = struct('node',node,'elem',elem);
option.L0 = 2;
option.maxIt = 4;
option.printlevel = 1;
option.plotflag = 0;
option.elemType = 'Q1';

%% Non-empty Dirichlet boundary condition.
pde = sincosdata;
mesh.bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
femPoisson(mesh,pde,option);

%% Pure Neumann boundary condition.
pde = sincosNeumanndata;
mesh.bdFlag = setboundary(node,elem,'Neumann');
femPoisson(mesh,pde,option);

%% Pure Robin boundary condition.
option.plotflag = 0;
pde = sincosRobindata;
mesh.bdFlag = setboundary(node,elem,'Robin');
femPoisson(mesh,pde,option);

%% Conclusion
% To do:
% 1. Fix the getH1errorQ1. The H1 error is computed wrong.
% 2. For pure Neuman boundary condition, the rate is not quite right.