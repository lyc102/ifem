%% CUBEAFEMQUADCURL quad curl equations on the unit cube
%
%   CUBEAFEMQUADCURL computes ND0-(CR-P0)-(linear ND) nonconforming approximations 
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
erruL2 = zeros(maxIt,1); 
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
[node,elem,bdFlag,HB] = uniformrefine3(node,elem,bdFlag,HB); 

%% PDE and options
pdeOrig = quadCurlDataSmooth1;
option.printlevel = 0;

%% AFEM cycles using only 1st Maxwell+ Stokes, 2nd Maxwell solved only for ref
for k = 1:maxIt
    fprintf('\n\n\nRefinement level %d\n', k);
    
    %% SOLVE the 1st Maxwell problem
    pdeCurl1.J = pdeOrig.quadcurlu;
    pdeCurl1.g_D = @(p) zeros(size(p,1),3);
    option.solver = 'mg';
    soln.w = Maxwellsaddle(node,elem,bdFlag,pdeCurl1,option,HB);
    
    %% SOLVE the Stokes problem
    
    % RHS of Stokes
    [curlwh, volume, curlbasis] = curlu3(node,elem,soln.w);
    pdeStokes.f = curlwh;
    pdeStokes.g_D = pdeOrig.curlu;
    option.solver = 'diag';
    option.printlevel = 1;
    solnStokes = Stokes3CRP0(node,elem,bdFlag,pdeStokes,option,HB);
    
    %% SOLVE the 2nd Maxwell problem (not used in the indicator for mesh refinement)
    
    % RHS os the 2nd Maxwell
    [elem2face,face] = dof3face(elem);
    NF = size(face,1); NT = size(elem,1); 
    soln.phi = reshape(solnStokes.u,NF,3);
    
    phihCenter = zeros(NT,3); % at each elem center
    for j = 1:3 % each component
        phihj = soln.phi(:,j);
        phihj2elem = phihj(elem2face);
        phihCenter(:,j) = sum(phihj2elem,2)/4;
        % Crouzeix-Raviart shape function 1-3\lambda = 1/4 at center
    end
    
    [elem2edge,edge] = dof3edge(sort(elem,2)); % sort for edge element
    NE = size(edge,1);
    bt = zeros(NT,6);
    for j = 1:6 % 6 basis of the edge element
        bt(:,j) = dot(curlbasis(:,:,j),phihCenter,2).*volume;
    end
    fCurl = zeros(2*NE,1);
    fCurl(1:NE,:) = accumarray(elem2edge(:),bt(:),[NE 1]);
    
    % solve the last Maxwell problem
%     option.solver = 'direct';
    option.solver = 'mg';
    pdeCurl2.J = fCurl;
    pdeCurl2.g_D = pdeOrig.exactu;
    pdeCurl2.g = pdeOrig.g;
    soln.u = Maxwell1saddle(node,elem,bdFlag,pdeCurl2,option,HB);
    %% ESTIMATE
    option.includeU = false; % only use the indicator by w_h and phi_h
    soln.curlw = curlwh;
    N(k) = 2*length(soln.u)+3*length(soln.phi);
    [etaK, est] = estimatequadcurl(node,elem,soln,pdeOrig,option);
    errphi(k) = getL2error3(node,elem,pdeOrig.curlu,soln.phi);
    etaTotal(k) = sqrt(sum(etaK.^2));
    curluh = curlu3(node,elem,soln.u);
    [erru(k), errK] = getL2error3(node,elem,pdeOrig.curlu,curluh);
    erruL2(k) = getL2error3ND1(node,elem,pdeOrig.exactu,soln.u);
    etaw(k) = sqrt(sum((est.eta1).^2));
    etaphi(k) = sqrt(sum((est.eta2).^2));
    etau(k) = sqrt(sum((est.eta3).^2));
    
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
    if N(k) > maxN; break; end
    
    markedElem = mark(elem,etaK,theta);
    if k < maxIt
        [node,elem,bdFlag,HB] = bisect3(node,elem,markedElem,bdFlag,HB);
    end
    
   
end

%%
close all;
figure(1);
r1 = showrate(N,erru,10,'k-o');
r2 = showrate(N,errphi,10,'r-o');
r3 = showrate(N,etaTotal,10,'b-*');
r4 = showrate(N,erruL2,10,'color',[0.8,0.2,0.7],'marker','*');
% r4 = showrate(N,etaw,10,'color',[0.2, 0.8,0.7],'marker','*');
% r5 = showrate(N,etaphi,10,'color',[0.2, 0.8,0.7],'marker','*');
% r6 = showrate(N,etau,10,'color',[0.2, 0.8,0.7],'marker','*');
set(gca,'xscale','log','yscale','log');
T = title('Convergence');
XL = xlabel('$\#$ DoFs');
YL = ylabel('Error Magnitudes');
L = legend('$\Vert  \nabla\times(u-u_h)\Vert$', ...
    ['$N^{' num2str(r1) '}$'],...
    '$\Vert  \nabla\times u - \phi_h\Vert$', ...
    ['$N^{' num2str(r2) '}$'],...
    '$\eta(w_h,\phi_h)$', ...
    ['$N^{' num2str(r3) '}$'],...
    '$\Vert u - u_h\Vert$', ...
    ['$N^{' num2str(r4) '}$'],'LOCATION','Best');
%     '$\eta_1(w_h)$', ['$N^{' num2str(r4) '}$'],...
%     '$\eta_2(\phi_h)$', ['$N^{' num2str(r5) '}$'],...
%     '$\eta_3(u_h)$', ['$N^{' num2str(r6) '}$'],...
    
set([XL,YL,L],'Interpreter','latex','FontSize', 16);
set(T,'Interpreter','latex','FontSize',20);
set(gca,'TickLabelInterpreter', 'latex');
grid on;
set(gcf,'color','w','Position', [100 200 500 800])
drawnow;

%% sanity check
if N<1e4
    figure(2);
    subplot(1,2,1);
    option.scale = 2;
    option.plot = 'quiver3';
    center = (node(elem(:,1),:) + node(elem(:,2),:) + ...
        node(elem(:,3),:) + node(elem(:,4),:))/4;
    showvector3(center,pdeOrig.curlu(center),option);
    subplot(1,2,2);
    curluh = curlu3(node,elem,soln.u);
    showvector3(center,curluh,option);
    
    set(gcf,'color','w','Position', [500 600 750 300])
    drawnow;
end