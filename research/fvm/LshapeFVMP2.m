function LshapeFVMP2
% Lsahpe solves Poisson equation in a L-shape domain with AFEM.
clear all; close all; clc;
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
theta = 0.5; % for bisection
maxIt = 50; % maxIt: maximum iteration number
N = zeros(maxIt,1); % store number of dof in each iteration 
err = zeros(maxIt,3); 
%%  Generate an initial mesh
node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
bdEdge = setboundary(node,elem,'Dirichlet');
elem = fixorientation(node,elem);   % counter-clockwise oritentation
elem = label(node,elem);   % label the mesh by the longest edge rule
%% pde data
pde.exactu = @exactu;
pde.f = @f;
pde.g_D = @g_D;
pde.Du = @Du;
pde.g_N=[];
pde.d=[];
%%  Get a fine mesh by uniform bisection
for k = 1:4
    [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
end
%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: Solve
    [u,A,~,~,edge] = PoissonP2FVM(node,elem,[],@f,@g_D);
    figure(1);  showresult(node,elem,u(1:size(node,1)),[-50,12]); 
    N(k) = size(u,1);
    % Step 2: Estimate
%    eta = estimaterecovery(node,elem,u(1:N(k))); % eta^2 when mark   
     eta = recoverflux2(node,elem,pde,bdEdge,u);
%    eta = recoverFlux(node,elem,pde,bdEdge);
    % Step 3: Mark
    markedElem = mark(elem,eta,theta);
%    markedElem = 1:size(elem,1);
    % compute error
   uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
   err(k,1) = sqrt((u-uI)'*A*(u-uI)); 
   err(k,2) = getL2error(node,elem,@exactu,u);
   err(k,3) = getH1error(node,elem,@pde.Du,u);
    % Step 4: Refine
    [node,elem] = bisect(node,elem,markedElem);
    bdEdge = setboundary(node,elem,'Dirichlet');
    if (N(k)>1e4), break; end   
end
N= N(1:k); 
erruIuh= err(1:k,1); errL2 = err(1:k,2); errH1 = err(1:k,3); 
figure(2);
r1 = showrate(N,erruIuh,20,'-*');
hold on;
r2 = showrate(N,errL2,20,'k-+');
r3 = showrate(N,errH1,20,'r-+');
legend('||Du_I-Du_h||',['N^{' num2str(r1) '}'], ...
       '||u_I-u_h||',['N^{' num2str(r2) '}'], ...
       '||Du-Du_h||',['N^{' num2str(r3) '}'], ...
       'LOCATION','Best');
return
%--------------------------------------------------------------------------
% Sub functions called by LSHAPE
%--------------------------------------------------------------------------
function z = f(p)
% load data (right hand side function)
z = zeros(size(p,1),1);
%--------------------------------------------------------------------------
function z = g_D(p) 
% Dirichlet boundary condition
z = exactu(p);
%--------------------------------------------------------------------------
function z = exactu(p) 
% exact solution of the test problem
r = sqrt(sum(p.^2,2));
theta = atan2(p(:,2),p(:,1));
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
z = r.^(2/3).*sin(2*theta/3);
%--------------------------------------------------------------------------
function z = u_x(p) 
% x-derivative of the exact solution
r = sqrt(sum(p.^2,2));
theta = atan2(p(:,2),p(:,1));
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
z = 2/3*r.^(-1/3).*cos(theta).*sin(2/3*theta) ...
  - 2/3*r.^(-1/3).*sin(theta).*cos(2/3*theta);
%--------------------------------------------------------------------------
function z = u_y(p) 
% y-derivative of the exact solution
r = sqrt(sum(p.^2,2));
theta = atan2(p(:,2),p(:,1));
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
z = 2/3*r.^(-1/3).*sin(theta).*sin(2/3*theta) ...
  + 2/3*r.^(-1/3).*cos(theta).*cos(2/3*theta);
function z = Du(p)
z = [u_x(p) u_y(p)];
%--------------------------------------------------------------------------