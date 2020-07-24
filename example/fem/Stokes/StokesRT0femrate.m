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
% In the error table, (u_h,p_h) is the RT0-P0 approximation of velocity and
% pressure. u_I is the canonical edge interpolant of RT0 element and p_I is
% the interpolant of p at the barycenter of each triangle. w = rotu is the
% vorticity, w_h = rot_hu_h is the numerical approximation and w_I is the
% Lagrange interpolation of w in P1 finite element space.
% 
% The rate of convergence of RT0-P0 depends crucially on the symmetry of
% the grids. For general unstructured grids, only half order for vorticity
% (i.e. velocity in the energy norm) and pressure. When the grid satisfies
% the approximately parallegram property, optimal first order is recovered.
% Furthermore for zero Dirichlet boundary condition, even second order for
% vorticity and pressure is observed.
%
% Simply change the method to BDM1b-P0 will lead to optimal order of
% convergence.
%
% See also
%  
%       StokesBDM1bfemrate    
%
% Created by Ming Wang, at Nov., 2011. Updated by Long Chen on Apr, 2013.

%% Example 1: Unit square with bisection grids
clear variables; 
close all;
% setup
[node,elem] = squaremesh([0,1,0,1],0.5);
pde = Stokesdata1; 
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% option
option.L0 = 1;
option.elemType = 'RT0-P0';
option.maxIt = 4;
option.solver = 'mg';
option.refType = 'bisect';
option.viewanglep = [-21,88];
option.viewanglew = [31,88];
% fem
femStokesHdiv(mesh,pde,option);
%% 
% For bisection grids, only half order for vorticity and pressure. It is
% interesting to note that the energy norm $||u_I - u_h||_1$ doesn't converge
% due to the loss of consistency $|| rot u - rot_h u_I||$ for the standard
% edge interpolant. Also no convergence of the maximum norm for pressure at the
% barycenter. It is observed vorticity oscillates on the boundary.

%% Example 2: Unit square with regular grids.
clear variables; 
close all;
% setup
[node,elem] = squaremesh([0,1,0,1],0.5);
pde = Stokesdata1;  % non-zero Dirichlet boundary condition
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% option
option.L0 = 1;
option.elemType = 'RT0-P0';
option.maxIt = 4;
option.solver = 'mg';
option.printlevel = 1;
option.refType = 'red';
option.viewanglep = [-21,88];
option.viewanglew = [31,88];
% fem
femStokesHdiv(mesh,pde,option);

%%
% The same pde data is used but the refinement rule is changed to the red
% refinement. For the red refinement grids and zero Dirichlet boundary
% condition, second order convergence for vorticity, and discrete pressure
% error (i.e. comparing to the interpolant p_I) is observed. The energy
% norm of u_I - u_h is 1.5 order.