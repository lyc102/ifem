close all; 
clear variables;
pde = HodgeLaplacian3Edata1;
[node,elem] = cubemesh([0,1,0,1,0,1],0.25);
% bdFlag = setboundary3(node,elem,'Neumann');
% Pure Neumann boundary condition doesn't work.
bdFlag = setboundary3(node,elem,'Dirichlet');
% bdFlag = setboundary3(node,elem,'Dirichlet','x==0','Neumann','~(x==0)');

err = zeros(4,1); N = zeros(4,1);
% option.solver = 'diag';
option.solver = 'tri';
option.printlevel = 2;
option.mg.Vit = 1;
option.mg.smoothingstep = 3;    % Smoothing step.
option.mg.smoothingratio = 1.5; % ratio of variable smoothing
% option.solver = 'direct';
for i = 1:4
    [sigma,u] = HodgeLaplacian3E(node,elem,bdFlag,pde,option);
    err(i) = getL2error3(node,elem,pde.sigma,sigma);
    N(i) = size(u,1);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end
showrate(N,err,2);