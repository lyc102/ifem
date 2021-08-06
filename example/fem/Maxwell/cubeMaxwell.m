%% CUBEMAXWELL solves Maxwell type equations in a cube using lowest order element.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

close all;

%% Defacult setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);
pde = Maxwelldata2;
% pde = planewavedata1;
% bdFlag = setboundary3(node,elem,'Neumann');
bdFlag = setboundary3(node,elem,'Dirichlet');
option.solver = 'amg';
% option.printlevel = 2;

%% Parameters
maxIt = 4; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
uIuhErr = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    % when epsilon is zero, a mixed formulation is needed
    if(abs(pde.epsilon)>1.0e-8)
        [u,edge,eqn] = Maxwell(node,elem,bdFlag,pde,option); 
    else
        [u,edge,eqn] = Maxwellsaddle(node,elem,bdFlag,pde,option); 
    end
    % compute error
    tic;
    energyErr(k) = getHcurlerror3ND(node,elem,pde.curlu,real(u));
    L2Err(k) = getL2error3ND(node,elem,pde.exactu,real(u));
    uI = edgeinterpolate(pde.exactu,node,edge);
    uIuhErr(k) = sqrt((real(u)-uI)'*(eqn.A)*(real(u)-uI));        
    time = toc;
    fprintf('Time to compute the error %4.2g s\n',time);
    N(k) = length(u);
    h(k) = 1./(size(node,1)^(1/3)-1);   
end

%% Plot convergence rates
figure;
showrateh3(h,energyErr,1,'k-+','|| curl (u-u_h) ||',...
           h,uIuhErr,1,'r-+','|| curl (u_I-u_h) ||',...
           h,L2Err,1,'b-+','|| u-u_h||');