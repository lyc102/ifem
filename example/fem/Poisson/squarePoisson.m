%% SQUAREPOISSON Poisson equation in a square domain.
%
%   squarePoisson computes linear approximations of the Poisson equation in
%   the unit square on a sequence of meshes obtained by uniform refinement.
% 
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed. The 2nd order convergent rate between two
% discrete functions ||DuI-Duh|| is known as superconvergence.
%
% See also Poisson, crack, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; 
clear variables;

%% Parameters 
maxIt = 4; 
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1);

%% Generate an initial mesh 
% [node,elem] = squaremesh([0 1 0 1], 0.25);
load fourholes.mat
% bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
% bdFlag = setboundary(node,elem,'Neumann');
for k = 1:0
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
end

%% Get the data of the pde
pde = sincosdata;
% pde = mixBCdata;

%% Finite Element Method        
for k = 1:maxIt
    % refine mesh
   [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
    % solve the equation
    [soln,eqn,info] = Poisson(node,elem,bdFlag,pde);
%     [soln,eqn,info] = PoissonP2(node,elem,bdFlag,pde,option);
    uh = soln.u;
    % record and plot
    N(k) = size(elem,1);
    h(k) = 1./(sqrt(size(node,1))-1);
    if N(k) < 2e3 % show mesh and solution for small size
        figure(1);  showresult(node,elem,uh);    
    end
    % compute error
    uI = pde.exactu(node); % nodal interpolation
    errL2(k) = getL2error(node,elem,pde.exactu,uh);
    errH1(k) = getH1error(node,elem,pde.Du,soln.Du);
    erruIuh(k) = sqrt((uh-uI)'*eqn.Lap*(uh-uI));
end

%% Plot convergence rates and display error table
figure(2);
showrateh3(h,errH1,1,'-*','||Du-Du_h||',...
           h,errL2,1,'k-+','||u-u_h||', ...
           h,erruIuh,1,'m-+','||DuI-Du_h||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||DuI-Du_h||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e',erruIuh,'%0.5e');