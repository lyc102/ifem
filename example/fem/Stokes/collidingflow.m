%% COLLIDINGFLOW colliding flow in a square domain
%
% Stokes equations on the square [0,1]^2. A simple model of colliding flow.
% The force f = 0, the velocity  u1 = 20xy^3, u2 = 5x^4-5y^4, p = 60x^2y -
% 20y^3.
%
% Dirichlet boundary condition is imposed.
%
% Reference: page 237 in Finite Elements and Fast Iterative Solvers with
% Applications in Incompressible Fluid Dynamics. by Howard C. Elman, David
% J. Silvester, and Andrew J. Wathen.
%
% See also squareStokes, femStokes
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; 
clear variables;

%% Setting
% mesh
[node,elem] = squaremesh([0 1 0 1], 0.25);
figure; showmesh(node,elem);
% pde
pde = Stokesdata1; 
bdFlag = setboundary(node,elem,'Dirichlet');
% options
option.L0 = 0;
option.maxIt = 4;

%% P2-P1 element
disp('P2-P1')
option.elemType = 'P2P1';
femStokes(node,elem,bdFlag,pde,option);

%% P2-P0 element
disp('P2-P0')
option.elemType = 'P2P0';
femStokes(node,elem,bdFlag,pde,option);

%% isoP2-P0 element
disp('isoP2-P0')
option.elemType = 'isoP2P0';
femStokes(node,elem,bdFlag,pde,option);

%% isoP2-P1 element
disp('isoP2-P1')
option.elemType = 'isoP2P1';
femStokes(node,elem,bdFlag,pde,option);

%% P1b-P1 element
disp('P1b-P0')
option.elemType = 'P1bP1';
femStokes(node,elem,bdFlag,pde,option);

%% CR-P0 element
disp('CR-P0')
option.elemType = 'CRP0';
femStokes(node,elem,bdFlag,pde,option);