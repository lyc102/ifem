%%
close all; 
clear variables;
[node,elem] = squaremesh([0,1,0,1],1/4);
pde = HodgeLaplacianEdata1;
% pde = HodgeLaplacianFdata1;

%% Edge element space
option.elemType = 'ND0';
% option.solver = 'direct';
option.solver = 'diag';
option.L0 = 2;
option.mg.Vit = 1;
option.mg.smoothingstep = 2;    % Smoothing step.
% option.mg.smoothingratio = 1.5; % ratio of variable smoothing
bdFlag = setboundary(node,elem,'Dirichlet','x==0','Neumann','~(x==0)');
mfemHodgeLap(node,elem,bdFlag,pde,option);