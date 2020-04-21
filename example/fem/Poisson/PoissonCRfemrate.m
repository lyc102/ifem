%% RATE OF CONVERGENCE OF NONCONFORMING LINEAR ELEMENT FOR POISSON EQUATION
%
% This example is to show the rate of convergence of CR non-conforming finite element
% approximation of the Poisson equation on the unit square with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
% - Robin boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem Poissonfemrate">Poissonfemrate</a>.
%
% See also PoissonP2femrate, Poissonafemrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
mesh = struct('node',node,'elem',elem);
option.L0 = 2;
option.maxIt = 4;
option.printlevel = 1;
option.elemType = 'CR';
option.plotflag = 1;

%% Non-empty Dirichlet boundary condition.
pde = sincosdata;
mesh.bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
femPoisson(mesh,pde,option);

%% Pure Neumann boundary condition.
option.plotflag = 0;
% pde = sincosNeumanndata;
pde = sincosdata;
mesh.bdFlag = setboundary(node,elem,'Neumann');
femPoisson(mesh,pde,option);

%% Pure Robin boundary condition.
pde = sincosRobindata;
mesh.bdFlag = setboundary(node,elem,'Robin');
femPoisson(mesh,pde,option);

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. No superconvergence for ||DuI-Duh||.
%
% MGCG converges uniformly in all cases.