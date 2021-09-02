[node,elem] = cubemesh([-1,1,-1,1,-1,1],1); 
bdFlag = setboundary(node,elem,'Dirichlet'); 
pde = mixBCdata3;
maxIt = 3;
option.solver = 'direct';

%% Solve and plot
err = zeros(maxIt,3); N = zeros(maxIt,1);
for i =1:maxIt
%     temperr = mixedPoisson(node,elem,pde,bdFlag,option);
    [u,sigma,eqn] = Poisson3RT0(node,elem,bdFlag,pde,option);
%     err(i,1) = temperr(1);
%     err(i,2) = temperr(2);
    err(i,1) = getL2error3RT0(node,elem,pde.Du,sigma);
    sigmaI = faceinterpolate3(pde.Du,node,elem);
    err(i,2) = getL2error3RT0(node,elem,pde.Du,sigmaI);
% %     err(i,3) = getL2error3RT0(node,elem,pde.Du,sigma,[]);
    err(i,3) = sqrt((sigma-sigmaI)'*eqn.M*(sigma-sigmaI));
%     err(i,4) = getHdiverror3RT0(node,elem,pde.f,-sigma,[]);
    N(i) = size(elem,1);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
end
figure(1); hold on; clf;
r3 = showrate3(N,err(:,1),1,'r-','L2sigma',N,err(:,2),1,'k-','L2u',...
               N,err(:,3),1,'b-','\sigma_I - \sigma_h');
