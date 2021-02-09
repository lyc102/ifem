%% LSHAPEAFEMQUADCURL quad curl equations on the Lshaped domain
%
%   LSHAPEAFEMQUADCURL computes ND0-(CR-P0)-(linear ND) nonconforming approximations 
%   of the quad curl equations in the unit cube. The mesh is refined
%   adaptively guided by a residual-based estimator using a separate marking
%   strategy with two marking parameters for eta_1 which estimates the error
%   for w and eta_2 for phi.
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; clear;

%% Set up
maxN = 6e5;
maxIt = 24;
theta = 0.3;
theta1 = 0.5;
nDofMaxwell = zeros(maxIt,1);
nDofStokes = zeros(maxIt,1);
h = zeros(maxIt,1);
erru = zeros(maxIt,1); 
erruL2 = zeros(maxIt,1); 
errw = zeros(maxIt,1);
errphiL2 = zeros(maxIt,1);
errphi = zeros(maxIt,1);
etaTotal = zeros(maxIt,1);
etaw = zeros(maxIt,1);
etaphi = zeros(maxIt,1);
etau = zeros(maxIt,1);


%% Generate initial mesh
load meshLshape3.mat % a minimal Lshape mesh
bdFlag = setboundary3(node,elem,'Dirichlet');
for i = 1
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end
[elem,bdFlag,HB] = label3(node,elem,'all',bdFlag);

%% PDE and options
gamma = 3+1/6; alpha = 2/3;
% gamma = 8/3; alpha = 2/3;
% gamma = 5/2; alpha = 1/2;
% gamma = 7/3; alpha = 1/3; % extreme case
pde = quadCurlDataLshape1(gamma,alpha);
% pde = quadCurlDataLshape2(gamma,alpha);
% pde = quadCurlDataLshape(gamma,alpha);
option.printlevel = 0;

%% AFEM cycles using only 1st Maxwell+ Stokes, 2nd Maxwell solved only for ref
for k = 1:maxIt
    fprintf('\n\n\nRefinement level %d\n', k);
    NT = size(elem,1); 
    
    %% SOLVE the 1st Maxwell problem
    
    pdeCurl1.J = pde.quadcurlu;
    pdeCurl1.g_D = @(p) zeros(size(p,1),3);
    option.solver = 'direct';
    if NT > 2e4; option.solver = 'diag'; end
    fprintf('\n**************Solving the first Maxwell problem for w**************');
    soln.w = Maxwellsaddle(node,elem,bdFlag,pdeCurl1,option,HB);
%     soln.w = Maxwellsaddle(node,elem,bdFlag,pdeCurl1,option);
    fprintf('*********************************************************************\n');
    %% SOLVE the Stokes problem
    
    % RHS of Stokes
    [curlwh, volume, curlbasis] = curlu3(node,elem,soln.w);
    pdeStokes.f = curlwh;
    pdeStokes.g_D = pde.curlu;
    option.solver = 'diag';
    option.printlevel = 1;
    fprintf('\n**************Solving the Stokes problem for phi*******************\n');
    solnStokes = Stokes3CRP0(node,elem,bdFlag,pdeStokes,option,HB);
    fprintf('*********************************************************************\n');
    %% SOLVE the 2nd Maxwell problem (not used in the indicator for mesh refinement)
    
    % RHS os the 2nd Maxwell
    [elem2face,face] = dof3face(elem);
    NF = size(face,1); 
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
    option.solver = 'direct';
    if NE > 2e4; option.solver = 'diag'; end
    pdeCurl2.J = fCurl;
    pdeCurl2.g_D = pde.exactu;
    pdeCurl2.g = pde.g;
    fprintf('\n*****Solving the second Maxwell problem for u (for reference)*****');
    soln.u = Maxwell1saddle(node,elem,bdFlag,pdeCurl2,option,HB);
    fprintf('*********************************************************************\n');
    %% ESTIMATE
    soln.curlw = curlwh;
    nDofMaxwell(k) = 2*NE;
    nDofStokes(k) = 3*NF;
    
    option.includeU = true; % only use the indicator by w_h and phi_h
    option.scaling = false; % no scaling and separate marking
    [etaK, est] = estimatequadcurl(node,elem,soln,pde,option);
    
    hK = est.hK;
%     h = sqrt(sum(hK.^2)/NT);
    Dlambda = gradbasis3(node,elem);
    curlphi = curlu3CR(elem2face,soln.phi,Dlambda);
    [errphi(k), errphiK] = getL2error3(node,elem,pde.curlcurlu,curlphi);
    [errphiL2(k), errKphiL2] = getL2error3(node,elem,pde.curlu,soln.phi);
    etaTotal(k) = sqrt(sum(etaK.^2));
    curluh = curlu3(node,elem,soln.u);
    [erru(k), errK] = getL2error3(node,elem,pde.curlu,curluh);
    erruL2(k) = getL2error3ND1(node,elem,pde.exactu,soln.u);
%     est.eta1 = est.eta1.*hK;
%     est.eta1 = est.eta1*min(hK);
    etaw(k) = sqrt(sum((est.eta1).^2));
    etaphi(k) = sqrt(sum((est.eta2).^2));
    etau(k) = sqrt(sum((est.eta3).^2));
    
    
    fprintf('Total eta          = %g \n', etaTotal(k));
%     fprintf('Mesh size          = %g \n', max(hK));
    fprintf('# face             = %g \n', NF);
    fprintf('# edge             = %g \n', NE);
    fprintf('||phi_h - curl u|| = %g \n', errphiL2(k));
    %% Visualize
    figure(1);
    subplot(1,3,1)
    etawE2V = accumarray(elem(:),repmat(est.eta1,[4,1]),[size(node,1), 1]);
    showsolution3(node,elem,etawE2V,'z<0.25 ','EdgeColor','k','facecolor','flat');
    colorbar;
    title('eta for w')
    view(120,45)
    
    subplot(1,3,2)
    
    % errphiL2E2V = accumarray(elem(:),repmat(errphiL2K,[4,1]),[size(node,1), 1]);
    % showsolution3(node,elem,errphiL2E2V,'z<0.25 ','EdgeColor','k','facecolor','flat');

    etaphiE2V = accumarray(elem(:),repmat(est.eta2,[4,1]),[size(node,1), 1]);
    showsolution3(node,elem,etaphiE2V,'z<0.25 ','EdgeColor','k','facecolor','flat');
    colorbar;
    title('eta for phi')
    view(120,45)
    
    subplot(1,3,3)
%     errU2v = accumarray(elem(:),repmat(erruK,[4,1]), [size(node,1), 1]);
%     showsolution3(node,elem,errU2v,'z<0.25','EdgeColor','k','facecolor','flat');
    % errL2U2v = accumarray(elem(:),repmat(abs(erruL2K),[4,1]),[size(node,1), 1]);
    % showsolution3(node,elem,errL2U2v,'z<0.25 ','EdgeColor','k','facecolor','flat');
    errphiE2V = accumarray(elem(:),repmat(errphiK,[4,1]),[size(node,1), 1]);
    showsolution3(node,elem,errphiE2V,'z<0.25 ','EdgeColor','k','facecolor','flat');
    title('H^1 error for phi')
    colorbar;
    view(120,45)
    
    set(gcf,'color','w','Position', [100 600 1400 300])
    drawnow;
    
    %% MARK and REFINE
    if nDofStokes(k) > maxN; break; end
%     
%     markedElem = mark(elem,etaK,theta);
%     markedElemC = mark(elem,est.elemRescurlwh,2/3*theta);
%     markedElem = unique([markedElem; markedElemC]);

    % markedElem = mark(elem,errKphiL2,theta); % cheating
    markedElem = mark2(elem,est.eta1,theta1,est.eta2,theta);
    
    fprintf('eta1(marked)^2/eta1^2 = %g \n', sum(est.eta1(markedElem).^2)/sum(est.eta1.^2));
    fprintf('eta2(marked)^2/eta2^2 = %g \n', sum(est.eta2(markedElem).^2)/sum(est.eta2.^2));
    if k < maxIt
        [node,elem,bdFlag,HB] = bisect3(node,elem,markedElem,bdFlag,HB);
    end
    
   
end

%%
close all;
figure(1);
nDofMaxwell = nDofMaxwell(1:k);
nDofStokes = nDofStokes(1:k);
erru = erru(1:k);
errphiL2 = errphiL2(1:k);
etaTotal = etaTotal(1:k);
erruL2 = erruL2(1:k);
etaw = etaw(1:k);
etaphi = etaphi(1:k);
etau = etau(1:k);
errphi = errphi(1:k);
ITER_THRESH = ceil(k/2);
r1 = showrate(nDofMaxwell,erru,ITER_THRESH,'k-o');
r2 = showrate(nDofStokes,errphiL2,ITER_THRESH,'r-o');
% r3 = showrate(nDofMaxwell+nDofStokes,etaTotal,ITER_THRESH,'b-*');
r3 = showrate(nDofMaxwell,erruL2,ITER_THRESH,'b-*');
r4 = showrate(nDofMaxwell/2,etaw,ITER_THRESH,'color',[0.7,0.2,0.7],'marker','*');
r5 = showrate(nDofStokes,etaphi,ITER_THRESH,'color',[0.2,0.8,0.7],'marker','*');
r6 = showrate(nDofMaxwell,etau,ITER_THRESH,'color',[0.7,0.8,0.2],'marker','*');
r7 = showrate(nDofStokes,errphi,ITER_THRESH,'g-o');
set(gca,'xscale','log','yscale','log');
T = title('Convergence');
XL = xlabel('$\#$ DoFs');
YL = ylabel('Error Magnitudes');
L = legend('$\Vert  \nabla\times(u-u_h)\Vert$', ...
    ['$(\# \mbox{edge})^{' num2str(r1) '}$'],...
    '$\Vert  \nabla\times u - \phi_h\Vert$', ...
    ['$(\# \mbox{face})^{' num2str(r2) '}$'],...
        '$\Vert u - u_h\Vert$', ...
    ['$(\# \mbox{edge})^{' num2str(r3) '}$'],...
     '$\eta_1(w_h)$', ['$(\# \mbox{edge})^{' num2str(r4) '}$'],...
    '$\eta_2(w_h, \phi_h)$', ['$(\# \mbox{face})^{' num2str(r5) '}$'],...
    '$\eta_3(\phi_h, u_h)$', ['$(\# \mbox{edge})^{' num2str(r6) '}$'],...
    '$|\nabla\times u-\phi_h|_{1,h}$', ['$(\# \mbox{face})^{' num2str(r7) '}$'],...
'LOCATION','southwest');
   
%     '$\eta(w_h,\phi_h)$', ...
%     ['$(\# \mbox{DoF})^{' num2str(r3) '}$'],...
    
set([XL,YL,L],'Interpreter','latex','FontSize', 16);
set(T,'Interpreter','latex','FontSize',20);
set(gca,'TickLabelInterpreter', 'latex');
grid on;
% set(gcf,'color','w','Position', [100 200 300 600])
set(gcf,'color','w','Position',[613 598 820 820])
drawnow;

%% sanity check
if nDofMaxwell<1e4
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
    
    set(gcf,'color','w','Position', [500 600 750 300])
    drawnow;
end