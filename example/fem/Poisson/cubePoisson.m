%% CUBEPOISSON solves Poisson equation in a cube using linear element.
%
% See also  
%   cubePoissonQ1, cubePoissonP2
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;
%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.25);

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
for k = 1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);  
%     [node,elem,HB,bdFlag] = uniformbisect3(node,elem,HB,bdFlag);  
    % solve the equation
    [soln,eqn] = Poisson3(node,elem,bdFlag,pde); 
    uh = soln.u;
    % compute error
    N(k) = size(node,1);
    h(k) = 1./(size(node,1)^(1/3)-1);    
    errH1(k) = getH1error3(node,elem,pde.Du,soln.Du);
    errL2(k) = getL2error3(node,elem,pde.exactu,uh);
    uI = pde.exactu(node);  % nodal interpolation
    erruIuh(k) = sqrt((uh-uI)'*eqn.Lap*(uh-uI));
end

%% Plot convergence rates
figure;
showrateh3(h,errH1,2,'-*', '$|| Du - Du_h ||', ...
           h,errL2,2,'k-+', '|| u - u_h ||', ...
           h,erruIuh,2,'m-+','|| D u_I - D u_h ||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|| u-u_h ||','|| Du-Du_h ||','|| Du_I-Du_h ||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e',erruIuh,'%0.5e');       

%% Conclusion
%
% The error in H1 and L2 norm converges at optimal rate. And the nodal
% interpolant and the FEM solution is super-close which is known as 
% superconvergence of uniform meshes.