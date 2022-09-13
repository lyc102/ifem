%% test curl-curl equation in 2d
clc; clear;
%%
[node,elem] = squaremesh([-1,1,-1,1],0.25);
pde = struct();
pde.f_x = @(p) (1 - 2*pi^2)*sin(pi*p(:,2)).*cos(pi*p(:,1));
pde.f_y = @(p) (2*pi^2 - 1)*sin(pi*p(:,1)).*cos(pi*p(:,2));
pde.exactu = @(p) [sin(pi*p(:,2)).*cos(pi*p(:,1)) -sin(pi*p(:,1)).*cos(pi*p(:,2))];
pde.curlu = @(p) -2*pi*cos(pi*p(:,1)).*cos(pi*p(:,2));
pde.g_D = @(p) zeros(size(p,1),1);

%%
bdFlag = setboundary(node,elem,'Dirichlet');
uh = CurlCurl2dNd0(node,elem,bdFlag,pde);