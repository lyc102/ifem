%% Setting of the problem
global s
pde = checkboarddata; % f = (x1^2-1)*(x2^2-1); Kellogg problem
pde.L = 1;
option.theta = 0.3;
option.estType = 'star';
option.maxIt = 25;
option.maxN = 2e4;
option.solver = 'mg';
option.tol = 1e-6;
[node,elem] = squaremesh([-1,1,-1,1],0.5);
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
save Kelloggcf err1 err2 err3 err4
load Kelloggcf
figure; 
plot_error_table(err1.N,err1.energyError,err2.N,err2.energyError,...
                 err3.N,err3.energyError,err4.N,err4.energyError);
%Save figure
saveas(gcf, 'error_Kellogg_cf', 'pdf') 
figure;
plot_error_table(err1.N,err1.eta,err2.N,err2.eta,err3.N,err3.eta,err4.N,err4.eta);
ylabel('Estimator','interpreter','latex', 'FontSize', 22)
saveas(gcf, 'eta_Kellogg_cf', 'pdf') 