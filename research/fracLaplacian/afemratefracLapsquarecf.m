close all; clear all;

%% Setting of the problem
global s
pde = fracLapdata2; % f =  sin(2*pi*x1)*sin(2*pi*x2);
pde.L = 1;
option.theta = 0.3;
option.estType = 'star';
option.maxIt = 18;
option.maxN = 2e4;
option.solver = 'mg';
% option.solver = 'direct';
option.tol = 1e-6;
[node,elem] = squaremesh([0 1 0 1],0.25);
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
% save squarecf err1 err2 err3 err4
load squarecf
figure; 
plot_error_table(err1.N,err1.H1,err2.N,err2.H1,err3.N,err3.H1,err4.N,err4.H1);
%Save figure
saveas(gcf, 'error_smooth', 'pdf') 
figure;
plot_error_table(err1.N,err1.eta,err2.N,err2.eta,err3.N,err3.eta,err4.N,err4.eta);
ylabel('Estimator','interpreter','latex', 'FontSize', 22)
saveas(gcf, 'eta_smooth', 'pdf') 