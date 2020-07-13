%% HODGELAPMGRATE test multigrid solver for Hodge Laplacian
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
close all; clear variables;
pde = HodgeLaplacianEdata1;
% bdFlag = setboundary(node,elem,'Dirichlet','x==0','Neumann','~(x==0)');
option.elemType = 'ND0';

%% Options
% option.N0 = 10;  % number of nodes in the coarsest mesh
option.maxIt = 4;
option.rateflag = 0;
% option.mg.solver = 'VCYCLE';    % V,W;
% option.mg.coarsematrix = 'G'; % Galerkin formulation
option.mg.smoothingstep = 2;    % Smoothing step.
option.mg.smoothingratio = 1.5; % ratio of variable smoothing
option.mg.Vit = 1;      % number of cycles for Schur complement   

%% Square
disp('Square Domain')
option.L0 = 3;  % refine to a finer mesh
[node,elem] = squaremesh([0,1,0,1],1/4);
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

%% Diagonal Preconditioner
option.solver = 'sdiag';
mfemHodgeLap(mesh,pde,option);

%% Triangular Preconditioner
option.solver = 'tri';
mfemHodgeLap(mesh,pde,option);

%% Approximated Factorization Preconditioner
option.solver = 'appf';
mfemHodgeLap(mesh,pde,option);

%% Lshape domain
disp('Lshape Domain')
option.L0 = 2;  % refine to a finer mesh
pde = fveconedata;
[node,elem] = squaremesh([-1,1,-1,1],0.25);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

%% Diagonal Preconditioner
option.solver = 'sdiag';
mfemHodgeLap(mesh,pde,option);

%% Triangular Preconditioner
option.solver = 'tri';
mfemHodgeLap(mesh,pde,option);

%% Approximated Factorization Preconditioner
option.solver = 'appf';
mfemHodgeLap(mesh,pde,option);


%% Crack domain
disp('Crack Domain')
pde = fveconedata;
node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];        % nodes
elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];            % elements
elem = label(node,elem);                        % label the mesh
[node,elem ] = uniformrefine(node,elem);
bdFlag = setboundary(node,elem,'Dirichlet');
option.L0 = 3;  % refine to a finer mesh
mesh = struct('node',node,'elem',elem,'bdFlag',bdFlag);

%% Diagonal Preconditioner
option.solver = 'sdiag';
mfemHodgeLap(mesh,pde,option);

%% Triangular Preconditioner
option.solver = 'tri';
mfemHodgeLap(mesh,pde,option);

%% Approximated Factorization Preconditioner
option.solver = 'appf';
mfemHodgeLap(mesh,pde,option);