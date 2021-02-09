%% CUBEAFEMQUADCURL quad curl equations on the unit cube
%
%   CUBEAFEMQUADCURL computes ND0-(CR-P0)-ND0 nonconforming approximations 
%   of the quad curl equations in the unit cube. The mesh is refined
%   adaptively guided by a residual-based estimator.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear;

%% Set up
maxN = 5e5;
maxIt = 20;
theta = 0.3;
N = zeros(maxIt,1);
h = zeros(maxIt,1);
erru = zeros(maxIt,1); 
errw = zeros(maxIt,1);
errphi = zeros(maxIt,1);
etaTotal = zeros(maxIt,1);
etaw = zeros(maxIt,1);
etaphi = zeros(maxIt,1);
etau = zeros(maxIt,1);

%% Generate initial mesh
H = pi;
[node,elem,HB] = cubemesh([0,H,0,H,0,H],H/2);
bdFlag = setboundary3(node,elem,'Dirichlet');

%% PDE and options
pde = quadCurlDataSmooth1;

option.solver = 'direct';
option.printlevel = 0;

%% Finite Element Method        
for k = 1:maxIt
    
    %% SOLVE
    [soln,eqn] = quadcurl3NC(node,elem,bdFlag,pde,option);
    N(k) = 2*length(soln.u)+3*length(soln.phi);
    %% ESTIMATE
    [erru(k), errK] = getHcurlerror3ND(node,elem,pde.curlu,soln.u);
    [etaK, est] = estimatequadcurl(node,elem,soln,pde,option);
    etaTotal(k) = sqrt(sum(etaK.^2));
    etaw(k) = sqrt(sum((est.eta1).^2));
    etaphi(k) = sqrt(sum((est.eta2).^2));
    etau(k) = sqrt(sum((est.eta3).^2));
    fprintf("Iter %2d, # Dof= %6d; \n||curl(u-u_h)|| = %8.6g,eta = %8.6g \n \n",...
           k, N(k), erru(k), etaTotal(k))
    %% Visualize
    figure(1);
    set(gcf, 'Position', [100 600 1000 300])
    subplot(1,3,1)
    showboundary3(node,elem,'z<pi/2','facealpha',0.5);
    
    subplot(1,3,2)
    errorE2V = accumarray(elem(:),repmat(errK,[4,1]),[size(node,1), 1]);
    showsolution3(node,elem,errorE2V,'z<pi/2 ','EdgeColor','k');
    
    subplot(1,3,3)
    etaE2V = accumarray(elem(:),repmat(etaK,[4,1]),[size(node,1), 1]);
    showsolution3(node,elem,etaE2V,'z<pi/2 ','EdgeColor','k');
    
    drawnow;
    %% MARK and REFINE
    if N(k) > maxN
        break;
    end
    
    markedElem = mark(elem,etaK,theta);
    if k < maxIt
        [node,elem,bdFlag,HB] = bisect3(node,elem,markedElem,bdFlag,HB);
        [elem,bdFlag] = sortelem3(elem,bdFlag);
    end
    
end

%%
close all;
figure(1);
r1 = showrate(N,erru,10,'r-o');
r2 = showrate(N,etaTotal,10,'b-*');
r3 = showrate(N,etaw,10,'color',[0.2, 0.8,0.7],'marker','*');
r4 = showrate(N,etaphi,10,'color',[0.2, 0.8,0.7],'marker','*');
r5 = showrate(N,etau,10,'color',[0.2, 0.8,0.7],'marker','*');
set(gca,'xscale','log','yscale','log');
T = title('Convergence');
XL = xlabel('$\#$ DoFs');
YL = ylabel('Error Magnitudes');
L = legend('$\Vert  \nabla\times(u-u_h)\Vert$', ...
    ['$N^{' num2str(r1) '}$'],...
    '$\eta(w_h,\phi_h,u_h)$', ...
    ['$N^{' num2str(r2) '}$'],...
    '$\eta_1(w_h)$', ['$N^{' num2str(r3) '}$'],...
    '$\eta_2(\phi_h)$', ['$N^{' num2str(r4) '}$'],...
    '$\eta_3(u_h)$', ['$N^{' num2str(r5) '}$'],...
    'LOCATION','Best');
set([XL,YL,L],'Interpreter','latex','FontSize', 12);
set(T,'Interpreter','latex','FontSize',16);
grid on;
set(gcf,'color','w','Position', [100 200 350 500])
drawnow;


%% sanity check
if N<5e3
    figure(2);
    subplot(1,2,1);
    option.scale = 2;
    option.plot = 'quiver3';
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
        node(elem(:,3),:) + node(elem(:,4),:))/4;
    showvector3(center,pde.curlu(center),option);
    subplot(1,2,2);
    curluh = curlu3(node,elem,soln.u);
    showvector3(center,curluh,option);
    % showvector3(center,curlwh_comp2,option);
    
    set(gcf,'color','w','Position', [500 600 750 300])
    drawnow;
end