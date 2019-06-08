function cubePoissonP2
%% CUBEPOISSONP2 solves Poisson equation in a cube using quadratic element.
%
% Copyright (C) 2008 Long Chen. See COPYRIGHT.txt for details.

close all;
%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);

%% Get the data of the pde
pde = sincosdata3;
% pde = polydata3;

%% Set up boundary condition
bdFlag = setboundary3(node,elem,'Dirichlet');
% bdFlag = setboundary3(node,elem,'Neumann');

%% Finite Element Method        
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1);
for k=1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);  
    % solve the equation
    [soln,eqn] = Poisson3P2(node,elem,bdFlag,pde); 
    uh = soln.u;
    N(k) = length(uh);
    h(k) = 1./(size(node,1)^(1/3)-1);    
    errH1(k) = getH1error3(node,elem,pde.Du,uh);
    errL2(k) = getL2error3(node,elem,pde.exactu,uh);    
    uI = pde.exactu([node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2]);
    erruIuh(k) = sqrt((uh-uI)'*eqn.Lap*(uh-uI));
    %error(k) = norm(u(1:N(k))-uI(1:N(k)));
end

%% Plot convergence rates
figure;
showrateh3(h,errH1,2,'-*', '|| Du - Du_h ||', ...
           h,errL2,2,'k-+', '|| u - u_h ||', ...
           h,erruIuh,2,'m-+','|| D u_I - D u_h ||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|| u-u_h ||','|| Du-Du_h ||','|| Du_I-Du_h ||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e',erruIuh,'%0.5e');       