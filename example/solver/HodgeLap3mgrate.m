%% HODGELAP3MGRATE test multigrid solver for Hodge Laplacian in 3D
%
% Reference 
%
% L. Chen, Y. Wu, L. Zhong and J. Zhou. Multigrid Preconditioners for Mixed
% Finite Element Methods of Vector Laplacian. Journal of Scientific
% Computing, 77:101?128, 2018. https://doi.org/10.1007/s10915-018-0697-7
%
% See also HodgeLap3mgrate
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

%%
close all;
clear option
option.elemType = 'ND0';

%% Options
option.maxIt = 4;
option.rateflag = 0;
% option.mg.solver = 'VCYCLE';    % V,W;
% option.mg.coarsematrix = 'G'; % Galerkin formulation
option.mg.smoothingstep = 3;    % Smoothing step.
option.mg.smoothingratio = 1.5; % ratio of variable smoothing
option.mg.Vit = 1;      % number of cycles for Schur complement 
option.mg.x0 = 'rand';

%% Cube
pde = HodgeLaplacian3Edata1;
[node,elem] = cubemesh([0,1,0,1,0,1],0.25);
bdFlag = setboundary3(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

%% Diagonal Preconditioner
option.solver = 'sdiag';
mfemHodgeLap3(mesh,pde,option);

%% Triangular Preconditioner
option.solver = 'tri';
mfemHodgeLap3(mesh,pde,option);

%% Approximated Factorization Preconditioner
option.solver = 'appf';
mfemHodgeLap3(mesh,pde,option);

%% Lshape domain
pde = fveconedata;
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);
% [node,elem] = delmesh(node,elem,'x>0 & y<0');
[node,elem] = delmesh(node,elem,'x<0 & y<0 & z>0');
showboundary3(node,elem);
bdFlag = setboundary3(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

%% Diagonal Preconditioner
option.solver = 'sdiag';
mfemHodgeLap3(mesh,pde,option);

%% Triangular Preconditioner
option.solver = 'tri';
mfemHodgeLap3(mesh,pde,option);

%% Approximated Factorization Preconditioner
option.solver = 'appf';
mfemHodgeLap3(mesh,pde,option);
