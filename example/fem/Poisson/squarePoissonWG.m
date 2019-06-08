%% SQUAREPOISSONWG Poisson equation in a square domain solved by Weak Galerkin method.
%
%   squarePoissoWGn computes the lowest order Weak Galerkin approximations
%   of the Poisson equation in the unit square on a sequence of meshes
%   obtained by uniform refinement. 
%
% See also PoissonWG, 
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear variables;

%% Parameters 
maxIt = 4; 
N = zeros(maxIt,1);
h = zeros(maxIt,1);
erruIuh = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem] = squaremesh([0 1 0 1], 0.25);
%bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
bdFlag = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');
%bdFlag = setboundary(node,elem,'Dirichlet');
for k = 1:3
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
end

%% Get the data of the pde
pde = sincosdata;
% pde = mixBCdata;
% option.solver = 'mg';

%% Finite Element Method        
for k = 1:maxIt
    % solve the equation
    [soln,eqn,info] = PoissonWG(node,elem,bdFlag,pde);
    N(k) = size(elem,1);
    h(k) = 1./(sqrt(size(node,1))-1);
    % compute error
    uI = zeros(N(k)+size(eqn.edge,1),1);
    uI(1:N(k)) = pde.exactu((node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3);
    uI(N(k)+1:end) = pde.exactu((node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2);
    erruIuh(k) = sqrt((soln.u-uI)'*eqn.A*(soln.u-uI));
    % refine mesh
   [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
end

%% Plot convergence rates
figure(2);
showrateh(h,erruIuh);