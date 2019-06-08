%% RATE OF CONVERGENCE OF MIXED FINITE ELEMENT METHOD (BDM1-P0) FOR POISSON EQUATION
%
% This example is to show the rate of convergence of RT0-P0 mixed finite
% element approximation of the Poisson equation on the unit square with the
% following boundary conditions:
%
% - pure Dirichlet boundary condition.
% - Pure Neumann boundary condition.
% - mxied boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem PoissonRT0femrate">PoissonRT0femrate</a>.
%
% The optimal rates of convergence for u and sigma are observed, namely,
% 1st order for L2 norm of u, H(div) norm of sigma, and 2nd order for 
% L2 norm of sigma. 
%
% The 2nd order convergent rates between two discrete functions ||uI-uh|| 
% is known as superconvergence.
%
% The 3rd order convergent rates between two discrete functions 
% ||sigmaI-sigmah|| for the Pure Neumman boundary condition 
% is known as superconvergence. The 2nd accurate numerical quadrature
% is required for the integral of rhs to observe such superconvergence.
%
% Triangular preconditioned GMRES (default solver)  and Uzawa
% preconditioned CG converges uniformly in all cases. Traingular
% preconditioner is two times faster than PCG although GMRES is used.
%
% See also PoissonRT0mfemrate, Poissonfemrate, Poissonafemrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
mesh = struct('node',node,'elem',elem);
option.L0 = 2;
option.maxIt = 4;
option.printlevel = 1;
option.elemType = 'BDM1';
option.solver = 'tri';
pde = sincosNeumanndata;

%% Pure Neumann boundary condition.
mesh.bdFlag = setboundary(node,elem,'Neumann');
mfemPoisson(mesh,pde,option);

%% Pure Dirichlet boundary condition.
mesh.bdFlag = setboundary(node,elem,'Dirichlet');
mfemPoisson(mesh,pde,option);

%% Mix Dirichlet and Neumann boundary condition.
option.solver = 'uzawapcg';
mesh.bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
mfemPoisson(mesh,pde,option);