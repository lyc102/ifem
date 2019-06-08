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
% option
option.L0 = 1;
option.elemType = 'BDM1B-P0';
option.maxIt = 4;
option.solver = 'mg';
option.refType = 'bisect';
% fem
femStokesHdiv(node,elem,bdFlag,pde,option);
%% 
% For bisection grids, only half order for vorticity and pressure. It is
% interesting to note that the energy norm $||u_I - u_h||_1$ doesn't converge
% due to the loss of consistency $|| rot u - rot_h u_I||$ for the standard
% edge interpolant. Also no convergence of maximum norm for pressure at the
% barycenter.

%% Example 2: Unit square with regular grids. Non-zero Dirichlet condition.
clear variables; 
close all;
% setup
[node,elem] = squaremesh([0,1,0,1],0.5);
pde = Stokesdata1;  % non-zero Dirichlet boundary condition
bdFlag = setboundary(node,elem,'Dirichlet');
% option
option.L0 = 1;
option.elemType = 'BDM1B-P0';
option.maxIt = 4;
option.solver = 'mg';
option.printlevel = 1;
option.refType = 'red';
% fem
femStokesHdiv(node,elem,bdFlag,pde,option);
%%
% The same pde data is used but the refinement rule is changed to red
% refinement. For red refinement grids and zero Dirichlet boundary
% condition, second order convergence for vorticity, and discrete pressure
% error (i.e. comparing to the interpolant p_I) is observed. The energy
% norm of u_I - u_h is 1.5 order.