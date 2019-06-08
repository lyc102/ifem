function KelloggFVMP2
%% KEOLLOGG Problem
%
% KELLOGG solves a diffusion equation with jump coefficients with AFEM.
%
% KELLOGG(maxN,theta) solves the problem within maxN number of vertices. The
% input argument theta is a parameter used in the marking step. 
%
% The KELLOGG command, if no input arguments, use maxN = 5e3 and theta = 0.5. 
%
% EXAMPLE
%
%    Kellogg 
%
% See also  crack, Lshape
%
% Created by Ming Wang
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;clc;clear all;
%% Problem setting
% $$-\nabla\cdot(d\nabla u) = 0  \quad \Omega=(-1,1)\times (-1,1)$$ 
%
% $$u = g_D \quad \partial \Omega$$
%
% The diffusion constant is discontinous; see the figure below. We set a2 =
% 1; a1 = 161.4476387975881 and choose boundary condition g_D such that the
% exact solution is $z = r^{0.1}\mu(\theta)$ in the poloar coordinate, where
% the formula of mu can be found in exactu function.

% [x,y] = meshgrid(-1:1:1,-1:1:1); 
% z = 0*x;
% surf(x,y,z,'linewidth',2); view(2);
% axis equal; axis tight;
% text(0.5,0.5,'a1','FontSize',12,'FontWeight','bold');
% text(-0.5,-0.5,'a1','FontSize',12,'FontWeight','bold');
% text(-0.5,0.5,'a2','FontSize',12,'FontWeight','bold');
% text(0.5,-0.5,'a2','FontSize',12,'FontWeight','bold');

%% Parameters
maxN = 3e3;     theta = 0.5;    maxIt = 200; 
N = zeros(maxIt,1); errH1 = zeros(maxIt,1); errL2 = zeros(maxIt,1);
erruIuh = zeros(maxIt,1);
%%  Generate an initial mesh
node = [-1 -1; 1 -1; 1 1; -1 1];
elem = [2 3 1; 4 1 3];
for i = 1:2
	[node,elem] = uniformbisect(node,elem);
end
bdEdge = setboundary(node,elem,'Dirichlet');
%% Set up PDE data
pde.f = @f;
pde.g_D = @exactu;
pde.exactu = @exactu;
pde.d = @d;
pde.Du=@Du;
%%  Adaptive Finite Element Method
% *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*

for k = 1:maxIt
    % Step 1: SOLVE
    [u,A,~,~,edge] = PoissonP2FVM(node,elem,[],pde.f,pde.g_D,[],[],pde.d);
   % [u,edge,A] = PoissonP2(node,elem,[],@pde.f,@pde.g_D,[],[],@pde.d);
   %  [u,edge,A] = PoissonP2FVMc(node,elem,[],pde.f,pde.g_D,[],[],pde.d);
    % Plot mesh and solution
    figure(1);  showresult(node,elem,u(1:size(node,1)),[27,26]);    
    % Step 2: ESTIMATE
   % eta = estimaterecovery(node,elem,u);            % recovery type
    eta = recoverflux2(node,elem,pde,bdEdge,u);   
    % Record error and number of vertices
   uI = pde.exactu([node; (node(edge(:,1),:)+node(edge(:,2),:))/2]);
   erruIuh(k) = sqrt((u-uI)'*A*(u-uI)); 
   errL2(k) = getL2error(node,elem,@pde.exactu,u);
   errH1(k) = getH1error(node,elem,@pde.Du,u,pde.d);
    N(k) = size(node,1);
    if (N(k)>maxN), break; end        
    % Step 3: MARK
    markedElem = mark(elem,eta,theta);
    % Step 4: REFINE
    [node,elem,bdEdge] = bisect(node,elem,markedElem,bdEdge);
end

%%  Plot convergent rates in energy norm
N= N(1:k); 
errL2 = errL2(1:k);
errH1 = errH1(1:k);
erruIuh = erruIuh(1:k);
figure(2);
r1 = showrate(N,errL2,100,'-*');
hold on;
r2 = showrate(N,errH1,100,'r-+');
r3 = showrate(N,erruIuh,100,'k-+');
legend('||u-u_h||',['N^{' num2str(r1) '}'], ...
       '||Du-Du_h||_K',['N^{' num2str(r2) '}'], ...
       '||Du_I-Du_h||_K',['N^{' num2str(r3) '}'], ...
       'LOCATION','Best');
   return;
   
% End of function KELLOGG
%%  Data of Kellogg Problem
function z = f(p)
% load data (right hand side function)
z = zeros(size(p,1),1);
function z = d(p)  % Diffusion constant
idx = (p(:,1).*p(:,2) >0);
z  = ones(size(p,1),1);
z(idx) = 161.4476387975881;
function z = exactu(p) % exact solution
gamma = 0.1;
sigma = -14.9225565104455152;
rho = pi/4;
r = sqrt(sum(p.^2,2));
theta = atan2(p(:,2),p(:,1));
theta = (theta>=0).*theta + (theta<0).*(theta+2*pi);
mu = (theta>=0 & theta<pi/2).*cos((pi/2-sigma)*gamma).*cos((theta-pi/2+rho)*gamma)...
    +(theta>=pi/2 & theta<pi).*cos(rho*gamma).*cos((theta-pi+sigma)*gamma)...
    +(theta>=pi & theta<1.5*pi).*cos(sigma*gamma).*cos((theta-pi-rho)*gamma)...
    +(theta>=1.5*pi & theta<2*pi).*cos((pi/2-rho)*gamma).*cos((theta-1.5*pi-sigma)*gamma);
z = r.^gamma.*mu;
%--------------------------------------------------------------------------
function z=Du(p)
% the gradient of exact solution
gamma=0.1;
sigma=-14.92256510455152;
rho=pi/4;
theta = atan2(p(:,2),p(:,1));%jiaodu
theta = (theta>=0).*theta +(theta<0).*(theta+2*pi);
t=1+(p(:,2).^2)./(p(:,1).^2);
r = sqrt(sum(p.^2,2));
rg=r.^gamma;

ux1=(p(:,1)>=0.0 & p(:,2)>=0.0).*(rg.*gamma./r.*cos((pi/2-sigma)*gamma)./r.*p(:,1)...
    .*cos((theta-pi/2+rho)*gamma)...
    +rg.*cos((pi/2-sigma)*gamma).*sin((theta-pi/2+rho)*gamma)...
    *gamma.*p(:,2)./(p(:,1).^2)./t);

uy1=(p(:,1)>=0.0 & p(:,2)>=0.0).*(rg*gamma./r.*cos((pi/2-sigma).*gamma)...
    .*cos((theta-pi/2+rho).*gamma)./r.*p(:,2)...
    -rg.*cos((pi/2-sigma).*gamma).*sin((theta-pi/2+rho)*gamma)...
    *gamma./p(:,1)./t);

ux2=(p(:,1)<=0.0 & p(:,2)>=0.0).*(r.^(-1.9).*p(:,1)*gamma.*cos(rho.*gamma)...
    .*cos((theta-pi+sigma).*gamma)...
    +rg.*cos(rho*gamma).*sin((theta-pi+sigma).*gamma)*gamma...
    .*p(:,2)./(p(:,1).^2)./t);
uy2=(p(:,1)<=0.0 & p(:,2)>=0.0).*( r.^(-1.9).*p(:,2)*gamma.*cos(rho.*gamma)...
    .*cos((theta-pi+sigma)*gamma)...
    -rg.*cos(rho*gamma).*sin((theta-pi+sigma).*gamma)*gamma./p(:,1)./t);

ux3=(p(:,1)<=0.0 & p(:,2)<=0.0).*(r.^(-1.9).*p(:,1).*gamma.*cos(sigma.*gamma)...
    .*cos((theta-pi-rho).*gamma) ...
    +rg.*cos(sigma.*gamma).*sin((theta-pi-rho).*gamma).*gamma...
    .*p(:,2)./(p(:,1).^2)./t);
uy3=(p(:,1)<=0.0 & p(:,2)<=0.0).*(r.^(-1.9).*p(:,2)*gamma.*cos(sigma*gamma)...
    .*cos((theta-pi-rho).*gamma) ...
    -rg.*cos(sigma*gamma).*sin((theta-pi-rho).*gamma)*gamma./p(:,1)./t);

ux4=(p(:,1)>=0.0& p(:,2)<=0.0).*(r.^(-1.9).*p(:,1).*gamma.*cos((pi/2-rho)*gamma)...
    .*cos((theta-3*pi/2-sigma)*gamma) ...
     +rg.*cos((pi/2-rho).*gamma).*sin((theta-3*pi/2-sigma)*gamma)...
     *gamma.*p(:,2)./(p(:,1).^2)./t);
 
 uy4=(p(:,1)>=0.0 & p(:,2)<=0.0).*(r.^(-1.9).*p(:,2)*gamma.*cos((pi/2-rho)*gamma)...
    .*cos((theta-3*pi/2-sigma)*gamma)...
    -rg.*cos((pi/2-rho)*gamma).*sin((theta-3*pi/2-sigma)*gamma)...
    *gamma./p(:,1)./t);     

z(:,1)= ux2+ux1+ux3+ux4;
z(:,2)= uy2+uy1+uy3+uy4;
