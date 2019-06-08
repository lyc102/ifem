function crackFVMP2
%% CRACK Problem
%
% crack solves the Poisson equation $-\Delta u =f$ in $\Omega$ and $u =
% g_D$ on $\partial \Omega$ in a crack domain
% $\Omega=\{|x|+|y|<1\}\backslash \{0<= x <=1, y=0\}$
%  using adaptive finite element method (AFEM). We choose f=1 and g_D such
%  that the exact solution is $u = r^{\beta}\sin(\beta\theta)-0.25r^2,
%  \beta = 1/2$ in the polar coordinate.
%
% EXAMPLE
%
%    crack 
%
% See also  crack_performance, Lshape
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all;clc;clear;
%% Parameters
maxN = 2e3;     theta = 0.5;    maxIt = 50; 
N = zeros(maxIt,1); err = zeros(maxIt,3);    
%%  Generate an initial mesh
node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];        % nodes
elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];            % elements
elem = label(node,elem);                        % label the mesh
bdEdge = setboundary(node,elem,'Dirichlet');    % Dirichlet boundary condition
showmesh(node,elem);                            % plot mesh
findelem(node,elem);                            % plot element indices
findnode(node,1:5);                             % plot node indices
findedge(node,[1 5; 5 6],1);                    % plot the crack edge
text(node(6,1)+0.05,node(6,2)+0.075,int2str(6), ...
     'FontSize',14,'FontWeight','bold');
%% 
% node 1 and node 6 are the same point (1,0)

%%  Get a fine mesh by uniform bisection
for k = 1:2
    [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
end
clf; showmesh(node,elem);
%% Set up PDE data
pde.f = @f;
pde.exactu = @exactu;
pde.g_D = @exactu;
pde.Du=@Du;
pde.d=[];
%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
   [u,A,~,~,edge] = PoissonP2FVM(node,elem,[],@pde.f,@pde.g_D);
    % Plot mesh and solution
    figure(1);  showresult(node,elem,u(1:size(node,1)),[-7,12]);    
    % Step 2: ESTIMATE
   % eta = estimaterecovery(node,elem,u); % recovery type % need eta^2 for mark
    eta = recoverFlux2(node,elem,pde,bdEdge,u);
    % Record error and number of vertices
   uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
   err(k,1) = sqrt((u-uI)'*A*(u-uI)); 
   err(k,2) = getL2error(node,elem,@exactu,u);
   err(k,3) = getH1error(node,elem,@pde.Du,u);
    N(k) = size(node,1);
    if (N(k)>maxN), break; end        
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
    [node,elem,bdEdge] = bisect(node,elem,markedElem,bdEdge);
end

%% Plot convergence rates
N= N(1:k); 
erruIuh= err(1:k,1); errL2 = err(1:k,2); errH1 = err(1:k,3); 
figure(2);
r1 = showrate(N,erruIuh,5,'-*');
hold on;
r2 = showrate(N,errL2,5,'k-+');
r3 = showrate(N,errH1,10,'r-+');
legend('||Du_I-Du_h||',['cN^{' num2str(r1) '}'], ...
       '||u-u_h||',['cN^{' num2str(r2) '}'], ...
       '||Du-Du_h||',['cN^{' num2str(r3) '}'], ...
       'LOCATION','Best');
%%
% Using AFEM, we obtain optimal rate of convergence the error in the energy
% norm ($N^{-0.5}$) and in the $L^2$ norm ($N^{-1}$).
end % End of function CRACK


%%  Data of CRACK
function z = f(p)   % load data (right hand side function)
z = ones(size(p,1),1);
end

function z = exactu(p) % exact solution
r = sqrt(sum(p.^2,2));
z = sqrt(0.5*(r-p(:,1)))-0.25*r.^2;
end

function z = Du(p) % derivative of the exact solution
r = sqrt(sum(p.^2,2));
z(:,1) = (p(:,1)./r-1)./sqrt(8*(r-p(:,1)))-0.5*p(:,1); % u_x
z(:,2) = p(:,2)./r./sqrt(8*(r-p(:,1)))-0.5*p(:,2);     % u_y
end