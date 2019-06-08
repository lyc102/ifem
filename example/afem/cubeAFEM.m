%% CUBEAFEM Problem
%
% CUBEAFEM solves the Poisson equation $-\Delta u =f$ in $\Omega$ and $u =
% g_D$ on $\partial \Omega$ in a cubic domain
% $\Omega=(-1,1)\times (-1,1) \times (-1,1)$
%  using adaptive finite element method (AFEM). We choose f=1 and g_D such
%  that the exact solution is $u = exp(-10*r^2)$ in the polar coordinate.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],0.5);
bdFlag = setboundary3(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag,'HB',HB);
pde = cubeafemdata;
option.L0 = 0;
option.maxIt = 10;
option.viewcut = '~(x>=0 & y>=0 & z>=0)'; 
option.viewangle = [122,18];
[err,time,solver,eqn] = afemPoisson3(mesh,pde,option);
%%
% Using AFEM, we obtain optimal convergent rate for the error in the energy
% norm and L2 norm.
% c = {'DOF','Error (H1)', 'Error (L2)'};
% makeHtmlTable([N errH1 errL2],[],[],c,[],[15 3 3]);