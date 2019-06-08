%% RATE OF CONVERGENCE OF ADAPTIVE FINITE ELEMENT METHOD USING WG ELEMENT
%
% This example is to show the rate of convergence of the lowest order weak
% Galerkin element approximation of the second order elliptic equation.
%
% Reference
%
% L. Chen, J. Wang and X. Ye. A posteriori error estimates for Weak
% Galerkin finite element methods for second order elliptic problems.
% Journal of Scientific Computing, 59(2), 496-511, 2014.

close all
clear variables
%% Lshape problem
% mesh
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary3(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% pde
pde = Lshapedata3;
% option
format shorte
option.L0 = 1;
option.maxIt = 50;
option.printlevel = 1;
option.plotflag = 1;
option.maxN = 2e4;
% AFEM
err = afemPoisson3(mesh,pde,option);
% plot
figure;
showrate2(err.N,err.H1,10,'b-*','|| Du-Du_h||',err.N,err.eta,10,'k-+','eta');
% latexerrtable(err.N,[err.H1 err.eta])

%% Jump coefficients
clear variables
% The diffusion coefficent is piecewise constant with large jump.
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | (x==-1)');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);
% pde
pde = jumpmgdata1;
global epsilon
epsilon = 1e-4;
% option
option.L0 = 3;
option.maxIt = 50;
option.printlevel = 1;
option.plotflag = 1;
option.maxN = 1e4;
option.viewcut = '~(x>=0 & y<=0 & z>=0)';
option.viewangle = [59,20];
% AFEM
err = afemPoisson3(mesh,pde,option);
% plot
figure;
showrate2(err.N,err.H1,10,'b-*','||Du-Du_h||',err.N,err.eta,10,'k-+','eta');
% Note no exact solution to compute the error. Code the energy later on.