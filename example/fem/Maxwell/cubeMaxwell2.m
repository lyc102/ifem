%% CUBEMAXWELL2 solves Maxwell type equations using quadratic element.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;

%% Defacult setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
pde = Maxwelldata2;
% pde = planewavedata1;
% bdFlag = setboundary3(node,elem,'Neumann');
bdFlag = setboundary3(node,elem,'Dirichlet');
% option.solver = 'cg';

%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
energyErr = zeros(maxIt,1);
uIuhErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    % solve the equation
    [u,T,eqn] = Maxwell2(node,elem,bdFlag,pde); 
    % compute error
    energyErr(k) = getHcurlerror3ND2(node,elem,pde.curlu,real(u));
    L2Err(k) = getL2error3ND2(node,elem,pde.exactu,real(u));
    uI = edgeinterpolate2(pde.exactu,node,T.edge,T.face,T.face2edge);
    uIuhErr(k) = sqrt((u-uI)'*eqn.A*(u-uI));
%     L2Err(k) = sqrt(abs((u-uI)'*eqn.M*(u-uI)));
%     energyErr(k) = sqrt((u-uI)'*eqn.A*(u-uI) + L2Err(k)^2);
    N(k) = length(u);
    h(k) = 1./(size(node,1)^(1/3)-1);            
end

%% Plot convergence rates
figure;
showrateh3(h,energyErr,1,'k-+','|| curl (u-u_h) ||',...
           h,uIuhErr,1,'r-+','|| curl (u_I-u_h) ||',...
           h,L2Err,1,'b-+','|| u-u_h||');