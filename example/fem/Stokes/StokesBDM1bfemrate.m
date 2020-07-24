%% RATE OF CONVERGENCE OF THE TMAC DISCRETIZATION FOR THE STOKES EQUANTIONS
%
% This example is to show the rate of convergence of RT0-P0 finite element
% approximation of the Stokes equations
%
% - grad div u + curl rot u + grad p  = f   in $\Omega$      
%
% - div u   = 0   in $\Omega$
%
% for the Dirichlet boundary condition:
%
% u $\cdot$ t   = g_t   on $\Gamma$   
%
% u $\cdot$ n   = g_n   on $\Gamma$   
%
% In the error table, (u_h,p_h) is the BDM1B-P0 approximation of velocity
% and pressure. u_I is the canonical edge interpolant of BDM1 element and
% p_I is the interpolant of p at the barycenter of each triangle. w = rotu
% is the vorticity, w_h = rot_hu_h is the numerical approximation and w_I
% is the Lagrange interpolation of w in P1 finite element space.
% 
% The rate of convergence of BDM1B-P0 is robust to the symmetry of the
% grids. For general unstructured grids, optimal first order convergence of
% velocity in H1 norm and for pressure in L2 norm is recovered. For
% vorticity  and for recovered velocity, the order is 1.5 for L2 norm.
%
% Reference
%
% L. Chen, M. Wang, and L. Zhong. Convergence Analysis of The Triangular
% MAC Scheme for Stokes Equations. Journal of Scientific Computing, 63(3),
% 716-744, 2015.

%% Example 1: Unit square with bisection grids
clear variables; 
close all;
% setup
[node,elem] = squaremesh([0,1,0,1],0.5);
pde = Stokesdata1; % zero Dirichlet boundary condition
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% option
option.L0 = 1;
option.elemType = 'BDM1B-P0';
option.maxIt = 4;
option.solver = 'mg';
option.printlevel = 1;
option.refType = 'bisect';
option.viewanglep = [-21,88];
option.viewanglew = [31,88];
% fem
femStokesHdiv(mesh,pde,option);
%% 

%% Example 2: Unit square with regular grids. Non-zero Dirichlet condition.
clear variables; 
close all;
% setup
[node,elem] = squaremesh([0,1,0,1],0.5);
pde = Stokesdata1; % zero Dirichlet boundary condition
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% option
option.L0 = 1;
option.elemType = 'BDM1B-P0';
option.maxIt = 4;
option.solver = 'mg';
option.printlevel = 1;
option.refType = 'red';
option.viewanglep = [-21,88];
option.viewanglew = [31,88];
% fem
femStokesHdiv(mesh,pde,option);