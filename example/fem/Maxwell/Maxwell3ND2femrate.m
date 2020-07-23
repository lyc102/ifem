%% RATE OF CONVERGENCE OF EDGE ELEMENT FOR MAXWELL EQUATIONS
%
% This example is to show the rate of convergence of quadratic edge finite
% element approximation of the Maxwell equation on the unit cube with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem Maxwell3ND2femrate">Maxwell3ND2femrate</a>.
%
% See also Maxwell3ND0femrate, Maxwell3ND2femrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear variables

%% Setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
mesh = struct('node',node,'elem',elem);
option.L0 = 0;
option.maxIt = 4;
option.elemType = 'ND2';
option.printlevel = 1;
option.plotflag = 1;

%% Dirichlet boundary condition.
fprintf('Dirichlet boundary conditions. \n');    
pde = Maxwelldata2;
bdFlag = setboundary3(node,elem,'Dirichlet');
femMaxwell3(mesh,pde,option);

%% Neumann boundary condition.
fprintf('Neumann boundary condition. \n');
option.plotflag = 0;
pde = Maxwelldata2;
mesh.bdFlag = setboundary3(node,elem,'Neumann');
femMaxwell3(mesh,pde,option);

%% Conclusion 
%
% Both the H(curl)-norm and the L2-norm is 2nd order. 

% MGCG using HX preconditioner converges uniformly in all cases.