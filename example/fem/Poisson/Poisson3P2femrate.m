%% RATE OF CONVERGENCE OF QUADRATIC ELEMENT FOR POISSON EQUATION
%
% This example is to show the rate of convergence of quadratic finite
% element approximation of the Poisson equation on the unit cube with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem Poisson3P2femrate">Poisson3P2femrate</a>.
%
% See also Poisson3femrate, Poisson3CRfemrate, Poissonfemrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear variables

%% Setting
[node,elem] = cubemesh([0,1,0,1,0,1],0.5); 
mesh = struct('node',node,'elem',elem);
% option.L0 = 1;
option.maxIt = 4;
option.elemType = 'P2';
option.printlevel = 1;
option.plotflag = 1;

%% Non-empty Dirichlet boundary condition.
fprintf('Mixed boundary conditions. \n');    
% pde = polydata3;
pde = sincosdata3;
mesh.bdFlag = setboundary3(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
femPoisson3(mesh,pde,option);

%% Pure Neumann boundary condition.
fprintf('Pure Neumann boundary condition. \n');
pde = sincosdata3Neumann;
option.plotflag = 0;
mesh.bdFlag = setboundary3(node,elem,'Neumann');
femPoisson3(mesh,pde,option);

% %% Robin boundary condition.
% fprintf('Robin boundary condition. \n');
% pde = sincosRobindata3;
% mesh.bdFlag = setboundary3(node,elem,'Robin');
% femPoisson3(mesh,pde,option);

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (2nd order) and L2-norm
% (3rd order) is observed. The 3rd order convergent rate between two
% discrete functions ||DuI-Duh|| is known as superconvergence.
%
% MGCG converges uniformly in all cases.
