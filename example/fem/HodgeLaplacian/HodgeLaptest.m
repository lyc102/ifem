close all; 
clear variables
pde = HodgeLaplacianEdata1;
% pde = HodgeLaplacianFdata1;
[node,elem] = squaremesh([0,1,0,1],1/8);
% bdFlag = setboundary(node,elem,'Neumann');
% Pure Neumann boundary condition doesn't work.
% bdFlag = setboundary(node,elem,'Dirichlet');
bdFlag = setboundary(node,elem,'Dirichlet','x==0','Neumann','~(x==0)');
err = zeros(4,1); N = zeros(4,1);
option.solver = 'direct';
for i = 1:4
    [sigma,u] = HodgeLaplacianE(node,elem,bdFlag,pde,option);
%     [sigma,u,AD] = HodgeLaplacianF(node,elem,pde,bdFlag,option);
    err(i) = getL2error(node,elem,pde.sigma,sigma);
    N(i) = size(u,1);
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
end
showrate(N,err,2);