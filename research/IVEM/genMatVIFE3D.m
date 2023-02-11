function matrix = genMatVIFE3D(pde,mesh,fem,meshI,femI)
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
S = globMatrixVIFE3DStiff(pde.A,mesh,fem,fem);
% This matrix generator is a little faster than doing it separately.
% Sx = globMatrixIFE3D(pde.A,mesh,femI,fem,[1,0,0],fem,[1,0,0]);
% Sy = globMatrixIFE3D(pde.A,mesh,femI,fem,[0,1,0],fem,[0,1,0]);
% Sz = globMatrixIFE3D(pde.A,mesh,femI,fem,[0,0,1],fem,[0,0,1]);
% S = Sx+Sy+Sz;

%% 2. Face Matrices

%% 3. Generate the Right Hand Side Vector
rhsF = globRHSVIFE3D(pde.f, mesh, fem, [0,0,0]);
rhs = rhsF + femI.b;
% rhsE = globRHSIFE3DFace(pde.A, pde.exactu, fem, femIF, 1);
% rhsP = globRHSIFE3DFace(pde.one, pde.exactu, fem, femIF, 0);

%% 4. Dirichlet Boundary Conditions
Atotal = S + femI.K + femI.S;

% if option.mass == 1 
%     M = globMatrixVIFE3DMass(pde.A,mesh,fem,fem);  
%     Atotal = Atotal + M +femI.M;
% end

node = meshI.node;
bcind = fem.bcind;
[bc,mapper] = boundaryNode3D(mesh,node,bcind);
tu = feval(pde.gD,node(:,1),node(:,2),node(:,3));
ub = tu;
ub(mapper) = 0;
rhsB = Atotal*ub;
A = Atotal(mapper,mapper);
f = rhs(mapper) - rhsB(mapper);

%% 5. Outputs
matrix = struct('A',A,'f',f,'S',S,'rhsF',rhsF,'tu',tu,'mapper',mapper);