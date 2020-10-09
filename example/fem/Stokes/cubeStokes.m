%% CUBESTOKES Stokes equations on the unit cube
%
%   SQUARESTOKE computes CR-P0 approximations of the Stokes equations in
%   the unit cube.
%
% Added by Shuhao Cao. Apr, 2020.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear;

%% Set up
maxIt = 4;
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
erruH1 = zeros(maxIt,1);
erruL2 = zeros(maxIt,1);
errp = zeros(maxIt,1);

%% Generate initial mesh
[node,elem] = cubemesh([-0.5,0.5,-0.5,0.5,-0.5,0.5],0.25);
bdFlag = setboundary3(node,elem,'Dirichlet');

%% PDE and options
pde.exactu = @(p)[(1-cos(2*pi*p(:,1))).*sin(2*pi*p(:,2)).*sin(pi*p(:,3)),... 
    -(1-cos(2*pi*p(:,2))).*sin(2*pi*p(:,1)).*sin(pi*p(:,3)), 0*p(:,3)];
pde.f = @(p) [-pi^2*(9*cos(2*pi*p(:,1))-5).*sin(2*pi*p(:,2)).*sin(pi*p(:,3))+ p(:,1).^2,... 
    pi^2*(9*cos(2*pi*p(:,2))-5).*sin(2*pi*p(:,1)).*sin(pi*p(:,3)), 0*p(:,3)];
pde.g = @(p) zeros(size(p,1),1);
pde.exactp = @(p) p(:,1).^3/3;
pde.g_D = @(p) pde.exactu(p);
% solver
option.solver = 'diag'; % diagonal preconditioner

%% Finite Element Method        
for k = 1:maxIt
    
    [soln,eqn] = Stokes3CRP0(node,elem,bdFlag,pde,option);
    uh = soln.u;
    uhvec = reshape(uh,size(uh,1)/3,3);
    ph = soln.p;
    N(k) = length(uh)+length(ph);
    h(k) = 1./((size(node,1)).^(1/3)-1);
    % compute error;
    uI = pde.exactu((node(eqn.face(:,1),:)+node(eqn.face(:,2),:)+node(eqn.face(:,3),:))/3);
    erruH1(k) = sqrt((uh-uI(:))'*eqn.A*(uh-uI(:)));
    erruL2(k) = getL2error3(node,elem,pde.exactu,uhvec);
    errp(k) = getL2error3(node,elem,pde.exactp,ph);
    fprintf('Number of Dof %d \n', N(k));
    
    if k < maxIt
        [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag); 
    end 
end

%% Plot convergence rates and display error
figure(2);
showrateh3(h,erruH1,1,'-*','| u_I-u_h |_1',...
           h,erruL2,1,'-*','|| u-u_h ||',...
           h,errp,1,'-+','|| p-p_h||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|u_I-u_h|_1','||u-u_h||','||p-p_h||'};
disptable(colname,N,[],h,'%0.3e',erruH1,'%0.5e',erruL2,'%0.5e',errp,'%0.5e');
