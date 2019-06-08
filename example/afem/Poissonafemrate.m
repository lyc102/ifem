%% RATE OF CONVERGENCE OF AFEM: P1 Linear Element for Poisson 
%
% This example is to show that: adaptive finite element methods can
% recovery the optimal convergent rate for elliptic equations even the
% solution is singular.
%
% # Lshape problem.
% # Kellogg problem.

% help Poissonafemrate;
close all;
clear variables

%% Lshape problem
[node,elem] = squaremesh([-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
pde = Lshapedata;
option.L0 = 1;
option.maxIt = 25;
option.maxN = 4e3;
option.printlevel = 1;
option.plotflag = 1;
option.viewangle = [-50,12];
format shorte

err = afemPoisson(mesh,pde,option);

figure;
showrate2(err.N,err.H1,10,'k-*','||Du-Du_h||',err.N,err.eta,10,'-+','eta');

%% Kellogg problem
[node,elem] = squaremesh([-1 1 -1 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
pde = Kelloggdata;
option.L0 = 1;
option.maxIt = 100;
option.maxN = 1e4;
option.theta = 0.2;
option.plotflag = 1;
option.rateshift = 30;
option.viewangle = [27,26];

err = afemPoisson(mesh,pde,option);

figure;
showrate2(err.N,err.H1,20,'k-*','||Du-Du_h||',err.N,err.eta,40,'-+','eta');
% latexerrtable(err.N,[err.H1 err.eta])