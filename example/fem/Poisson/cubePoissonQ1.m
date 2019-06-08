%% CUBEPOISSONQ1 Poisson equation in a CUBE domain using trilinear cube element 
%
%   cubePoissonQ1 computes trilinear finite element approximations of the
%   Poisson equation in the unit cube on a sequence of hex meshes obtained by
%   uniform refinement. It plots the approximation err vs the number of
%   degree of freedoms.
% 
% See also  
%   cubePoisson, cubePoissonP2
% 
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all;

%% Parameters 
maxIt = 4; 
N = zeros(maxIt,1);
h = zeros(maxIt,1);

%% Domain and pde
h0 = 1.0/2.0;
cube = [0 1 0 1 0 1];

pde = sincosdata3;
option.solver = 'amg';
% option.isoparametric = true;

%% Finite Element Method        
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1);  
erruIuh = zeros(maxIt,1);
for k = 1:maxIt
    [node,elem] = cubehexmesh(cube,h0/(2.0^k));
    bdFlag = setboundary3(node,elem,'Dirichlet','abs(z)>eps','Neumann','abs(z)<eps');
%     bdFlag = setboundary3(node,elem,'Dirichlet');
%     [uh,Du,eqn,info] = Poisson3Q1(node,elem,pde,bdFlag,option);
    [uh,Du,eqn,info] = Poisson3T1(node,elem,pde,bdFlag,option);
    N(k) = size(node,1);
    h(k) = 1./(size(node,1)^(1/3)-1);    
    uI = pde.exactu(node);
    erruIuh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
    errH1(k) = getH1error3Q1(node,elem,pde.Du,uh);  
    errL2(k) = getL2error3Q1(node,elem,pde.exactu,uh);        
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
