function LshapeFVM
%% LSHAPE Problem
%
% LSHAPE solves the Poisson equation $-\Delta u =f$ in $\Omega$ and $u =
% g_D$ on $\partial \Omega$ in a crack domain $\Omega=(-1,1)^2\backslash
% (0,1)\times (-1,0)$
%  using adaptive finite element method (AFEM). We choose f and g_D such
%  that the exact solution is $u = r^{\beta}\sin(\beta\theta), \beta = 2/3$
%  in the polar coordinate.
%
% To see the improvement using AFEM, run Lshape_uniform to get the
% convergent rate using uniform refinement.
%
% EXAMPLE
%    Lshape 
%
% See also  crack, Lshape_uniform
%
% Created by Chen-Song Zhang. Modified by Long Chen.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; 
%% Parameters
maxN = 3e3;     theta = 0.5;    maxIt = 50; 
N = zeros(maxIt,1);   energy = zeros(maxIt,1);  uIuhErrH1 = zeros(maxIt,1);
errH1 = zeros(maxIt,1);
%%  Generate an initial mesh
node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0]; % nodes
elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
bdEdge = setboundary(node,elem,'Dirichlet');
elem = fixorientation(node,elem);   % counter-clockwise oritentation
elem = label(node,elem);            % label the mesh by the longest edge rule
showmesh(node,elem);                % plot mesh
findelem(node,elem);                % plot element indices
findnode(node);                     % plot node indices

%%  Get a fine mesh by uniform bisection
for k = 1:1
    [node,elem,bdEdge] = uniformbisect(node,elem,bdEdge);
end

clf; showmesh(node,elem);
%% Set up PDE data
pde.f = @f;
pde.g_D = @exactu;
pde.Du=@Du; % used for recoverFlux;
pde.d=[];
%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
for k = 1:maxIt
    % Step 1: SOLVE
    [u,~,~,A] = PoissonFVM(node,elem,pde);
    % Plot mesh and solution
   figure(1);  showresult(node,elem,u,[-50,12]);    
    % Step 2: ESTIMATE
%    eta = estimaterecovery(node,elem,u);         % recovery type
%     eta = estimateresidual(node,elem,u,pde);    % residual type
%    eta = recoverflux(node,elem,pde,bdEdge,u);
      eta = effrecflux(node,elem,pde,bdEdge,u);   % local L2 project.
    errH1(k) = getH1error(node,elem,@Du,u);
    % Record error and number of vertices
    energy(k) = u'*A*u;
    uI = exactu(node);
    uIuhErrH1(k) = sqrt((uI-u)'*A*(uI-u));
    N(k) = size(node,1);
    if (N(k)>maxN), break; end        
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
    [node,elem,bdEdge] = bisect(node,elem,markedElem,bdEdge);
end
%% Plot convergence rates
N= N(1:k); 
uIuhErrH1 = uIuhErrH1(1:k);
errH1 = errH1(1:k);
energyError = sqrt(energy(1:k)-energy(k));
figure(2);
r1 = showrate(N,uIuhErrH1,3,'-*');
hold on;
r2 = showrate(N(1:k-1),energyError(1:k-1),3,'k-+');
r3 = showrate(N,errH1,3,'g-+');
legend('||Du_I-Du_h||',['N^{' num2str(r1) '}'],...
       'sqrt{E(u_k)-E(u_i)}',['N^{' num2str(r2) '}'],...
       '||Du-Du_h||',['N^{' num2str(r3) '}'],...
       'LOCATION','Best')
   %%
% In this example, since f=0, the Dirichlet energy of u is $\|u\|_A^2$. By
% the minimization of the Galerkin projection, we compute $\|u-u_i\|_A^2
% \approx \|u_k - u_i\|_A^2 = \|u_k\|_A^2 -\|u_i\|_A^2$.
%
% We also compute the energy norm between the nodal interpolation and the
% finite element approximation. It is shown that $\|u_I-u_h\|_A$ admits
% convergent rate more than optimal one $N^{-1/2}$. This is known as
% superconvergence. For a finite element function v, the squared energy
% norm can be computed as $\|v\|_A^2 = v'*A*v$.
% End of function LSHAPE

function z = f(p)
% load data (right hand side function)
z = zeros(size(p,1),1);
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
