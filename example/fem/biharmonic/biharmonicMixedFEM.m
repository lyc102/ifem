%% BIHARMONICMIXEDFEM solves the biharmonic equation
%
%     laplace^2 u = f;   [0,1]^2
%               u = g_D;  
%            Dn u = g_N.  
%
% by writing in a mixed form
%     
%       w = Lap u
%   Lap w = f
%
% The boundary condition u = g_D is imposed into the space (Dirichlet) and
% the boundary condition Dn u = g_N is imposed as a Newmann boundary
% condition. 
%
% The rate of convergence for u is optimal but for w is sub-optimal. For
% linear element, optimal order for w is also observed.
%
% Created by Jie Zhou. Clean up by Long Chen. Further clean up on
% biharmonic functions are needed and multigrid solvers should be included.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all
clear variables;

%% Parameters
maxIt = 4; 
N = zeros(maxIt,1);    
h = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);     
erruH1 = zeros(maxIt,1);
errwL2 = zeros(maxIt,1);     
errwH1 = zeros(maxIt,1);

%%  Generate an initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.25);
bdFlag = setboundary(node,elem,'Dirichlet');    % Dirichlet boundary condition

for k = 1:1
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
%     [node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
end

%% Set up PDE data
pde = biharmonicdata;

%% Finite Element Method        
for k = 1:maxIt
    % refine mesh
   [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);     
   [w,u] = biharmonicP1(node,elem,bdFlag,pde);
%    [w,u] = biharmonicP2(node,elem,bdFlag,pde);
%    [w,u] = biharmonicP3(node,elem,bdFlag,pde);
   N(k) = size(w,1)+size(u,1);
   h(k) = 1./(sqrt(size(node,1))-1);
   erruL2(k) = getL2error(node,elem,pde.exactu,u);
   erruH1(k) = getH1error(node,elem,pde.Du,u);
   errwL2(k) = getL2error(node,elem,pde.exactw,w);
   errwH1(k) = getH1error(node,elem,pde.Dw,w);
end
 
%% Plot convergence rates and display error table
figure;
subplot(1,2,1);
showrateh2(h,erruH1,1,'-*','|| Du - Du_h||',...
           h,erruL2,1,'k-+','|| u - u_h||');
subplot(1,2,2)
showrateh2(h,errwH1,1,'c-*','|| Dw - Dw_h||',...
           h,errwL2,1,'m-+','|| w - w_h||');       

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||Dw-Dw_h||','||w - w_h||'};
disptable(colname,N,[],h,'%0.3e',erruL2,'%0.5e',erruH1,'%0.5e',errwH1,'%0.5e',errwL2,'%0.5e');

