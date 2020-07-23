%% RATE OF CONVERGENCE OF EDGE ELEMENT FOR MAXWELL EQUATIONS
%
% This example is to show the rate of convergence of lowest order edge finite element
% approximation of the Maxwell equation on the unit cube with the
% following boundary conditions:
%
% - Non-empty Dirichlet boundary condition.
% - Pure Neumann boundary condition.
%
% The basis, data structure and numerical test is summarized in <a
% href="matlab:ifem Maxwell3ND0femrate">Maxwell3ND0femrate</a>.
%
% See also Maxwell3ND1femrate, Maxwell3ND2femrate
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear variables

%% Setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
mesh = struct('node',node,'elem',elem);
option.L0 = 1;
option.maxIt = 4;
option.elemType = 'ND0';
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
% The optimal rate of convergence of both the H(curl)-norm and L2-norm
% (1st order) is observed. The 2nd order convergent rate between two
% discrete functions ||DuI-Duh||_{\infty} is known as superconvergence.
%
% MGCG using HX preconditioner converges uniformly in all cases.
