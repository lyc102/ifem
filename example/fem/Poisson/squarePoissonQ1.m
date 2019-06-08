%% SQUAREPOISSONQ1 Poisson equation in a square domain using bilinear quad element 
%
%   squarePoissonQ1 computes bilinear finite element approximations of the
%   Poisson equation in the unit square on a sequence of quad meshes obtained by
%   uniform refinement. It plots the approximation err vs the number of
%   degree of freedoms.
% 
% See also squarePoisson
%
% Author: Huayi Wei < huayiwei1984@gmail.com>
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all;
clear variables;

%% Parameters 
maxIt = 4; 
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1);
erruIuh = zeros(maxIt,1);

%% Generate an initial mesh 
node = [0, 0; 1, 0; 1, 1; 0,1];
% node = [0, 0; 2 -1; 2.5, 0; 2, 1];
elem = [1,2,3,4];
for k = 1:3
    [node,elem] = uniformrefinequad(node,elem);
end

%% Get the data of the pde
pde = sincosdata;

%% Finite Element Method     
for k = 1:4
    % refine mesh    
    [node,elem] = uniformrefinequad(node,elem);
    % solve the equation    
    [soln,eqn,info] = PoissonQ1(node,elem,[],pde);
    uh = soln.u;
    % record and plot    
    N(k) = size(elem,1);
    h(k) = 1./(sqrt(size(node,1))-1);    
    if N(k) < 2e3 % show mesh and solution for small size
        figure(1);  showsolution(node,elem,uh);    
    end
    % compute error
    uI = pde.exactu(node); % nodal interpolation
    errH1(k) = getH1errorQ1(node,elem,pde.Du,uh);  
    errL2(k) = getL2errorQ1(node,elem,pde.exactu,uh);
    erruIuh(k) = sqrt((uh-uI)'*eqn.Lap*(uh-uI));
end

%% Plot convergence rates and display error table
figure(2);
showrateh3(h,errH1,1,'-*','||Du-Du_h||',...
           h,errL2,1,'k-+','||u-u_h||', ...
           h,erruIuh,1,'m-+','||DuI-Du_h||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||DuI-Du_h||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e',erruIuh,'%0.5e');