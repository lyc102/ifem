%% adaptive test for interface problem on a square

clear; close all
%% Kellogg
pde = simpleInterfaceData(1,5);
bdExpr = '(x == -1) | ( x==1 ) | ( y==-1 ) | (y == 1)';
[node,elem] = squarequadmesh([-1,1,-1,1],1/4);
%% 
maxN = 1e4;
maxIt = 50;
theta = 0.3;
tol = 0.05;
errQDuDuhTotal = zeros(maxIt,1);
errH1uhuITotal = zeros(maxIt,1);
errQuadTotal = zeros(maxIt,1);
etaTotal = zeros(maxIt,1);
stabTotal = zeros(maxIt,1);
etaStabTotal = zeros(maxIt,1);
N = zeros(maxIt,1);
option = [];
uH1seminorm = 1;

%% 
elem = num2cell(elem,2);
elemRefine = elem;
elem = fixorientationpoly(node,elem);
elemQuad = cell2mat(cellfun(@(x) fliplr(x(1:4)), elemRefine,'UniformOutput',false));
hangNode = [];
%% 
for k = 1:maxIt
    %% solve
    N(k) = size(node,1);
    elemVertexNumber = cellfun('length',elem);
    fprintf("\n# Nodes             = %5d\n", N(k));
    fprintf("max vertex         = %5d\n", max(elemVertexNumber));
    fprintf("# of hanging nodes = %5d\n", length(hangNode));
    N(k) = N(k) - length(hangNode);
    [u,A,assembleTime,solverTime] = PoissonVEM(node,elem,pde);    
    fprintf("Assembling+solving:  %10.5f s\n", assembleTime+solverTime);
    uI = pde.exactu(node);
    %% plot
    if size(elem,1) <= 1e4 % only plot meshes of small size
        figure(1);clf;
        subplot(1,3,1); showmeshpoly(node,elem);
        subplot(1,3,2); showsolutionpoly(node,elem,u); view(120,56);
        subplot(1,3,3); showsolutionpoly(node,elem,uI-u); view(120,56);
        axis([min(node(:,1))-0.05 max(node(:,1)) min(node(:,2))-0.05 max(node(:,2))])
        set(gcf,'color',[0.8, 0.8, 0.8],'Position', [100 800 1000 400])
        view(120,50); % Kellogg
        drawnow;
        pause(0.1);
    end
    errH1uhuITotal(k) = sqrt((u-uI)'*A*(u-uI)); % get the error in H1 norm
    
    
    %% estimate
    tic;
    [QDuh,area,~,~,stabK] = graduVEM(node,elem,u);
    center = (node(elemQuad(:,1),:)+node(elemQuad(:,2),:)+...
        node(elemQuad(:,3),:)+node(elemQuad(:,4),:))/4;
    QDu = pde.Du(center);
    errQDuDuh = sqrt(sum((QDu - QDuh).^2.*area,2));
    errQDuDuhTotal(k) = norm(errQDuDuh);
    errQuadTotal(k) = getH1errorQ1(node,elemQuad,pde.Du,u,pde.d(center));
    option.method = 'RT';
    option.d = pde.d(center);
    [etaR, ~, etaStabK] = estimaterecoveryVEM(node,elem,u,option);
    etaRes = estimateresidualVEM(node,elem,u,pde);
    etaK = sqrt(etaR.^2 + etaStabK.^2);
    etaTotal(k) = norm(etaK);
    stabTotal(k) = norm(stabK);
    etaStabTotal(k) = norm(etaStabK);
    estimateTime = toc;
    fprintf("Estimating        :  %10.5f s\n", estimateTime);
    %% MARK and REFINE
    if N(k) > maxN || errQuadTotal(k)/uH1seminorm<tol; break; end
    
    tic;
    markedElem = mark(elem,etaK,theta);
%     markedElem = mark(elem,etaR,theta);
    markedElem = unique([markedElem; find(elemVertexNumber>16)]);
    markTime = toc;
    fprintf("Marking           :  %10.5f s\n", markTime);
    
    if k < maxIt
        [node,elem, hangNode, elemRefine, timeInfo] = ...
            quadsectpoly(node,elem,markedElem, hangNode, elemRefine, bdExpr);
        elemQuad = cell2mat(cellfun(@(x) fliplr(x(1:4)), elemRefine,'UniformOutput',0));
        fprintf("Add new nodes     :  %10.5f s\n", timeInfo.addNewNodes);
        fprintf("Add new elements  :  %10.5f s\n", timeInfo.addNewElem);
        fprintf("Find hanging nodes:  %10.5f s\n", timeInfo.findHangNodes);
        fprintf("Add hang type 1   :  %10.5f s\n", timeInfo.addHangNodes1);
        fprintf("Add hang type 2   :  %10.5f s\n", timeInfo.addHangNodes2);
        fprintf("Remove dupe nodes :  %10.5f s\n", timeInfo.removeDupe);
%         [flag, ixErrElem] = checkpoly(node,elem);
%         if flag; break; end
    end
    
    %%
    fprintf("|u - u_{Q1}}|_1 = %10.5f\n", errQuadTotal(k));
    fprintf("||Q(Du - Duh)|| = %10.5f\n", norm(errQDuDuh));
    fprintf("eta             = %10.5f\n", etaTotal(k));
    fprintf("uh stab         = %10.5f\n", stabTotal(k));
    fprintf("flux stab       = %10.5f\n", etaStabTotal(k));
end

%% result
close all;
figure(1);
subplot(1,3,1)
showmeshpoly(node,elem);
subplot(1,3,2)
showsolutionpoly(node,elem,u);
view(120,50);
subplot(1,3,3)
showsolutionpoly(node,elem,uI);
view(120,50); % Kellogg

set(gcf,'color',[0.8, 0.8, 0.8],'Position', [100 600 1200 400])
drawnow;
%%
figure(2);
N = N(1:k);
rCut = 10; % rate cut off
errQDuDuhTotal = errQDuDuhTotal(1:k);
errH1uhuITotal = errH1uhuITotal(1:k);
errQuadTotal = errQuadTotal(1:k);
etaTotal = etaTotal(1:k);
stabTotal = stabTotal(1:k);
etaStabTotal = etaStabTotal(1:k);
r1 = showrate(N,errQDuDuhTotal,rCut,'r-o');
r3 = showrate(N,etaTotal,rCut,'color',[0.2, 0.8, 0.7],'marker','square');
r4 = showrate(N,stabTotal,rCut,'color',[0.7, 0.8, 0.2],'marker','diamond');
r5 = showrate(N,etaStabTotal,rCut,'color',[0.7, 0.2, 0.7],'marker','pentagram');
r6 = showrate(N,errQuadTotal,rCut,'color',[0.5, 0.2, 0.5],'marker','hexagram');
set(gca,'xscale','log','yscale','log');
T = title('Convergence');
XL = xlabel('$\#$ DoFs');
YL = ylabel('Error Magnitudes');
L = legend('$\Vert  \Pi\nabla(u-u_h)\Vert$', ['$N^{' num2str(r1) '}$'],...
    '$\eta(\Pi u_h)$', ['$N^{' num2str(r3) '}$'],...
    '$\Vert(I-\Pi) u_h \Vert_{\ell^2}$', ['$N^{' num2str(r4) '}$'],...
    '$\Vert(I-\Pi) \sigma_h \Vert_{\ell^2}$', ['$N^{' num2str(r5) '}$'],...
    '$|\!|\!| u-u_h|\!|\!|$', ['$N^{' num2str(r6) '}$'],...
    'LOCATION','Bestoutside');
set([XL,YL,L],'Interpreter','latex','FontSize', 16);
set(T,'Interpreter','latex','FontSize',16);
grid on;
set(gcf,'color','w','Position', [100 200 600 500])
drawnow;
