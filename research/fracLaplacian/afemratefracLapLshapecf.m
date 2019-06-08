close all; clear all;

%% Setting of the problem
global s
pde = fracLapdata9; % f = sin(2pi*x1)*sin(pi*x2);
option.theta = 0.3;
option.estType = 'star';
option.maxIt = 18;
option.maxN = 2e4;
option.solver = 'mg';
option.tol = 1e-6;
[node,elem] = squaremesh([-1,1,-1,1],0.5);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');

%% s = 0.2
s = 0.2; %#ok<*NASGU>
err1 = afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.4
s = 0.4;
err2 = afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.6
s = 0.6;
err3 = afemfracLap(node,elem,pde,bdFlag,option);

%% s = 0.8
s = 0.8;
err4 = afemfracLap(node,elem,pde,bdFlag,option);

%% save data and plot the table
save Lshapecf err1 err2 err3 err4
% load Lshapecf
% figure; 
% plot_error_table(err1.N,err1.energyError,err2.N,err2.energyError,...
%                  err3.N,err3.energyError,err4.N,err4.energyError);
% saveas(gcf, 'error_Lshapecf', 'pdf') 
% figure;
% plot_error_table(err1.N,err1.eta,err2.N,err2.eta,err3.N,err3.eta,err4.N,err4.eta);
% ylabel('Estimator','interpreter','latex', 'FontSize', 22)
% saveas(gcf, 'eta_Lshapecf', 'pdf') 