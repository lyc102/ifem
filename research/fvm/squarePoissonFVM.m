function [N,errL2,errH1,erruIuh] = squarePoissonFVM
%% SQUAREPOISSON Poisson equation in a square domain.
%
%   squarePoisson computes approximations of the Poisson equation in the
%   unit square on a sequence of meshes obtained by uniform refinement. It
%   plots the approximation error (in L2 norm or H1 norm) vs the number of
%   nodes.
% 
% See also Poisson, crack, Lshape

% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all;clc; clear;
%% Parameters 
maxIt = 7; N = zeros(maxIt,1);
errL2 = zeros(maxIt,1); errH1 = zeros(maxIt,1); erruIuh = zeros(maxIt,1);

%% Generate an initial mesh 
node = [0,0; 1,0; 1,1; 0,1];    % nodes
elem = [2,3,1; 4,1,3];          % elements
bdEdge = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
%  bdEdge = setboundary(node,elem,'Dirichlet');
for i=1:2
 [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
end

%% Get the data of the pde
 %pde = jumpdata1; % haven't test Neumann problem for this data.
%  pde = sincosdata;
pde = mixBCdata;
%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    [u,A] = PoissonFVM(node,elem,pde,bdEdge,'longedgecenter','direct');
%    figure(1);  showresult(node,elem,u);    
    N(k) = size(elem,1);
    % compute error
    uI = pde.exactu(node); % nodal interpolation
    errL2(k) = getL2error(node,elem,pde.exactu,u);
    errH1(k) = getH1error(node,elem,pde.Du,u);
    erruIuh(k) = sqrt((u-uI)'*A*(u-uI));
    [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
end

%% Plot convergence rates
figure(2);
r1 = showrate(N,errH1,3,'-*');
hold on;
r2 = showrate(N,errL2,3,'k-+');
r3 = showrate(N(2:end),erruIuh(2:end),2,'m-+');
legend('||Du-Du_h||',['N^{' num2str(r1) '}'], ...
       '||u-u_h||',['N^{' num2str(r2) '}'], ...
       '||Du_I-Du_h||',['N^{' num2str(r3) '}'], 'LOCATION','Best');
%% Error table
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. The 2nd order convergent rate between two
% discrete functions ||DuI-Duh|| is known as superconvergence.
