%% MAXWELLSADDLE3MGRATE test multigrid solvers for Maxwell equation in saddle point form
%
% Reference 
%
% L. Chen, Y. Wu, L. Zhong and J. Zhou. Multigrid Preconditioners for Mixed
% Finite Element Methods of Vector Laplacian. Journal of Scientific
% Computing, 77:101?128, 2018. https://doi.org/10.1007/s10915-018-0697-7
%
% See also HodgeLapmgrate, HodgeLap3mgrate
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

close all; 
clear variables
%% Problem
pde = Maxwellsaddledata;

%% Mesh
% cube
% [node,elem] = cubemesh([0,1,0,1,0,1],0.25);
% % bdFlag = setboundary3(node,elem,'Neumann'); % Pure Neumann boundary condition doesn't work.
% bdFlag = setboundary3(node,elem,'Dirichlet');
% Lshape
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
% [node,elem] = delmesh(node,elem,'x>0 & y<0');
[node,elem] = delmesh(node,elem,'x<0 & y<0 & z>0');
% showboundary3(node,elem);
bdFlag = setboundary3(node,elem,'Dirichlet');

%% Options
maxIt = 4;
% option.solver = 'diag';
% option.solver = 'tri';
option.solver = 'mg';
% option.solver = 'direct';
option.printlevel = 2;
option.mg.Vit = 1;
option.mg.smoothingstep = 3;    % Smoothing step.
option.mg.smoothingratio = 1.5; % ratio of variable smoothing
option.mg.x0 = 'rand';

%% Finite Element Methods
err = zeros(maxIt,1); 
h = zeros(maxIt,1);
for k = 1:maxIt
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    u = Maxwellsaddle(node,elem,bdFlag,pde,option);
    err(k) = getHcurlerror3NE(node,elem,pde.curlu,real(u));
    h(k) = 1./(size(node,1)^(1/3)-1);   
end
showrateh(h,err,2);