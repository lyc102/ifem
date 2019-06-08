function N = squarePoissonP2FVM
%% SQUAREPOISSONP2FVM solves Poisson equation in a square domain with 
% quadratic FVM.

%--------------------------------------------------------------------------
% Copyright Long Chen. Updated on Sep-28-2008
%--------------------------------------------------------------------------

close all; clear all;clc;
%% Parameters 
maxIt = 5; N = zeros(maxIt,1); errH1 = zeros(maxIt,1);

%% Generate initial mesh
node = [0,0; 1,0; 1,1; 0,1];    % nodes
elem = [2,3,1; 4,1,3];          % elements
% showmesh(node,elem);
% [tempvar,tempvar,edge]=dofRT0(elem);
% findedge(node,edge);
%findelem(node,elem);
 for i = 1:3
 	%[node,elem] = uniformrefine(node,elem);    
     [node,elem] = uniformbisect(node,elem);    
  end
%bdEdge = setboundary(node,elem,'Dirichlet');
%bdEdge=[];
bdEdge = setboundary(node,elem,'Dirichlet','all','Neumann','y==1');

%% pde data
pde=mixBCdata;
pde.g_R=[];

%% Finite volume Method 
for k = 1:maxIt
    u = PoissonP2FVM(node,elem,bdEdge,pde,'longedgecenter','mg');
    errH1(k) = getH1error(node,elem,pde.Du,u);
    N(k) = length(u);
    [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
    % [node,elem,bdEdge] = uniformrefine(node,elem,bdEdge);
end

%% Plot convergence rates
figure(2);
r1 = showrate(N,errH1,2,'-*');
legend('||Du-Du_h||',['N^{' num2str(r1) '}'], ...
        'LOCATION','Best');
end