function matrix = genMatVIFE3DmgHcurl(A,mesh,fem,femI)
%% Generate global matrices and load vector of PPIFEM for 3D ellipic eq
%     -div(A grad u)  = f,    x\in \Omega
%      where A is a piecewise constant on Omega^+ and Omega^-.
% INPUTS:
% pde --- given data function from equation, e.g. 
%         pde.A --- diffusion coefficient
%         pde.f --- right hand side function
%         pde.gD --- Dirichlet boundary value function 
%         pde.one --- constant function 1.
% mesh --- mesh structure. 
% fem --- global degree of freedom of FEM 
% femI --- quadrature info on interface cells
% femIF --- quadrature info on interface faces, required in PPIFE.
% PPtype --- partial penalty type possible value 
%            'N' : Nonsymmetric PPIFE (e = 1)
%            'S' : Symmetric PPIFE (e = -1)
%            'I' : Incomplete PPIFE (e = 0)
% OUTPUTS:
% matrix.S --- stiffness matrix (w/o boundary condition)
% matrix.E --- consistence matrix (w/o boundary condition)
% matrix.P --- penalty matrix (w/o boundary condition)
% matrix.A --- final FEM matrix (after boundary condition)
% matrix.rhsF --- load vector (w/o boundary condition)
% matrix.rhsE --- consistence vector (=0 if no interface face on boundary)
% matrix.rhsP --- penalty vector (=0 if no interface face on boundary)
% matrix.f --- final RHS matrix (after boundary condition)

% Last Modified: 08/07/2020 by Xu Zhang

%% 1. Stiffness Matrix
S = globMatrixVIFE3DStiff(A,mesh,femI,fem,fem);

Ndof = size(femI.p,1);
bdidx = zeros(Ndof,1); 
isBdEdge = true(Ndof,1);
isBdEdge(femI.mapper) = false;
bdidx(isBdEdge) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
A = T*S*T + Tbd;

%% 5. Outputs
matrix = struct('A',A);