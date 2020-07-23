%% RATE OF CONVERGENCE OF EDGE ELEMENT FOR MAXWELL EQUATIONS
%
% This example is to show the rate of convergence of linear edge finite element
% approximation of the Maxwell equation on the unit cube with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem Maxwell3ND1femrate">Maxwell3ND1femrate</a>.
%
% See also Maxwell3ND0femrate, Maxwell3ND2femrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear variables

%% Setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
mesh = struct('node',node,'elem',elem);
option.L0 = 1;
option.maxIt = 4;
option.elemType = 'ND1';
option.printlevel = 1;
option.plotflag = 1;

%% Dirichlet boundary condition.
fprintf('Dirichlet boundary conditions. \n');    
pde = Maxwelldata2;
bdFlag = setboundary3(node,elem,'Dirichlet');
femMaxwell3(mesh,pde,option);

%% Pure Neumann boundary condition.
fprintf('Neumann boundary condition. \n');
option.plotflag = 0;
pde = Maxwelldata2;
mesh.bdFlag = setboundary3(node,elem,'Neumann');
femMaxwell3(mesh,pde,option);

%% Conclusion 
%
% The H(curl)-norm is still 1st order but the L2-norm is improved to 2nd
% order. 

% MGCG using HX preconditioner converges uniformly in all cases.
