
h0 = 0.5;
[node,elem] = cubemesh([0,1,0,1,0,1],h0);
% bdFlag = setboundary3(node,elem,'Dirichlet');
bdFlag = setboundary3(node,elem,'Neumann');
% pde = mixBCdata3;
% pde = mixedPossiondata;
pde = sincosdata3;
maxIt = 3;
option.solver = 'uzawapcg';
% option.solver = 'tri';

%% Solve and plot
err = zeros(maxIt,2); 
h = zeros(maxIt,1);
h(1) = h0;
for i = 1:maxIt
%     err(i,1:2) = mixedPoisson(node,elem,pde,bdFlag,option);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    [u,sigma,eqn] = Poisson3RT0(node,elem,bdFlag,pde,option);
    err(i,1) = getL2error3RT0(node,elem,pde.Du,sigma);
    err(i,2) = getL2error3(node,elem,pde.exactu,u);
    sigmaI = faceinterpolate3(pde.Du,node,elem);
%     err(i,3) = getL2error3RT0(node,elem,pde.Du,sigmaI);
    err(i,3) = sqrt((sigma-sigmaI)'*eqn.M*(sigma-sigmaI));
%     err(i,4) = getHdiverror3RT0(node,elem,pde.f,-sigma,[]);
    h(i) = h0/2^i;
end
figure(1); hold on; clf;
showrateh3(h,err(:,1),1,'r-','|| \sigma-\sigma_h||',...
           h,err(:,2),1,'k-','|| u - u_h||',...
           h,err(:,3),1,'b-','|| \sigma_I - \sigma_h||');