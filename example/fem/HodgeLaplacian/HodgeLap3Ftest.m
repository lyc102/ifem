close all; clear all
pde = HodgeLaplacian3Fdata1;
[node,elem] = cubemesh([0,1,0,1,0,1],0.5);
% bdFlag = setboundary3(node,elem,'Neumann');
% Pure Neumann boundary condition doesn't work.
% bdFlag = setboundary3(node,elem,'Dirichlet');
bdFlag = setboundary3(node,elem,'Dirichlet','x==0','Neumann','~(x==0)');

err = zeros(3,1); N = zeros(3,1);
option = [];
for i = 1:3
    [sigma,u] = HodgeLaplacian3F(node,elem,pde,bdFlag,option);
    err(i) = getL2error3NE(node,elem,pde.exactsigma,sigma);
    N(i) = size(u,1);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end
showrate(N,err,2);