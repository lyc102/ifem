
%%
close all; 
clear variables;
[node,elem] = cubemesh([0,1,0,1,0,1],0.5);
pde = HodgeLaplacian3Edata1;
bdFlag = setboundary3(node,elem,'Dirichlet'); 
% bdFlag = setboundary(node,elem,'Dirichlet','x==0','Neumann','~(x==0)');

%% Edge element space
option.elemType = 'ND0';
% option.solver = 'direct';
option.solver = 'tri';
% option.L0 = 2;
mfemHodgeLap3(node,elem,bdFlag,pde,option);